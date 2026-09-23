import numpy as np, json, time, sys
from sklearn.kernel_ridge import KernelRidge
from numpy.linalg import svd
t0=time.time()
task = sys.argv[1]
man = json.load(open('results/r1_task_freeze_manifest.json'))
rng = np.random.default_rng(20260922)
if task == 'replogle':
    sp7 = '../SP-007-perturb-shadows-r0'
    k_idx = json.load(open(f'{sp7}/results/accum/index.json'))
    r_idx = json.load(open(f'{sp7}/results/accum_rpe1/index.json'))
    K = np.load('results/resp_k562.npy'); R = np.load('results/resp_rpe1.npy')
    k_genes = [l.strip() for l in open(f'{sp7}/results/gene_names.txt')]
    r_genes = json.load(open('/tmp/sp7_rpe1_genes.json'))
    genes = man['replogle']['shared_genes_cpm1']
    kc = np.array([k_genes.index(g) for g in genes]); rc = np.array([r_genes.index(g) for g in genes])
    kt = {t:i for i,t in enumerate(k_idx['primary'])}; rt = {t:i for i,t in enumerate(r_idx['strict'])}
    targets = man['replogle']['targets_shared_strict']
    Ks = np.stack([K[kt[t]][kc] for t in targets]); Rs = np.stack([R[rt[t]][rc] for t in targets])
    folds = man['replogle']['folds']
    directions = {'K562->RPE1': (Ks, Rs), 'RPE1->K562': (Rs, Ks)}
else:
    F = np.load('results/resp_frangieh.npy')
    f_idx = json.load(open('results/accum/index.json'))
    f_genes = json.load(open('/tmp/sp8_frangieh_genes.json'))
    genes = man['frangieh']['shared_genes_cpm1_all_conditions']
    gc = np.array([f_genes.index(g) for g in genes])
    targets = f_idx['targets']; tidx = {t:i for i,t in enumerate(f_idx['targets'])}
    Fsel = np.stack([F[tidx[t]][:, gc] for t in targets])   # (229, 3, G)
    conds = f_idx['conditions']
    folds = man['frangieh']['folds']
    directions = {}
    for d in range(3):
        src = [c for c in range(3) if c != d]
        X = np.concatenate([Fsel[:, c, :] for c in src], axis=1)
        directions[f'{conds[src[0]]}+{conds[src[1]]}->{conds[d]}'] = (X, Fsel[:, d, :])
def predict_models(Xtr, Ytr, Xte):
    n = Xtr.shape[0]
    preds = {}
    preds['context_mean'] = np.tile(Ytr.mean(0), (Xte.shape[0], 1))
    G_src = Xtr.shape[1] // (Xtr.shape[1] // Ytr.shape[1]) if Xtr.shape[1] % Ytr.shape[1] == 0 else None
    preds['xcontext_mean'] = Xte[:, :Ytr.shape[1]] if Xte.shape[1] >= Ytr.shape[1] else Xte
    # nested alpha via inner 4-fold on training targets
    def inner_score(model_fn, params):
        idx = rng.permutation(n); scores = {p: [] for p in params}
        for f in range(4):
            te = idx[f::4]; tr = np.setdiff1d(idx, te)
            for p in params:
                m = model_fn(p); m.fit(Xtr[tr], Ytr[tr])
                pr = m.predict(Xtr[te])
                scores[p].append(np.mean([np.corrcoef(pr[i], Ytr[te][i])[0,1] for i in range(len(te))]))
        return max(params, key=lambda p: np.mean(scores[p]))
    best_a = inner_score(lambda a: KernelRidge(kernel='linear', alpha=a), [0.1, 1.0, 10.0, 100.0])
    m = KernelRidge(kernel='linear', alpha=best_a).fit(Xtr, Ytr); preds['ridge'] = m.predict(Xte)
    best_g = inner_score(lambda g: KernelRidge(kernel='rbf', alpha=1.0, gamma=g), [1e-6, 1e-5, 1e-4])
    m = KernelRidge(kernel='rbf', alpha=1.0, gamma=best_g).fit(Xtr, Ytr); preds['kernel_ridge'] = m.predict(Xte)
    def lowrank(k):
        U, S, Vt = svd(Ytr - Ytr.mean(0), full_matrices=False)
        Uk, Sk, Vk = U[:, :k], S[:k], Vt[:k, :]
        W = KernelRidge(kernel='linear', alpha=1.0).fit(Xtr, Uk*Sk)
        return W, Vk, Ytr.mean(0)
    best_k = inner_score(lambda k: _LR(None, None, k), [5, 10, 20])
    W, Vk, mu = lowrank(best_k)
    preds['lowrank'] = W.predict(Xte) @ Vk + mu
    best_knn = inner_score(lambda k: _KNN(k), [3, 5, 10])
    preds['knn'] = _knn_pred(Xtr, Ytr, Xte, best_knn)
    return preds
class _LR:
    def __init__(self, X, Y, k): self.k=k; self.X=X; self.Y=Y
    def fit(self, x, y):
        U, S, Vt = svd(y - y.mean(0), full_matrices=False)
        self.Vk = Vt[:self.k]; self.mu = y.mean(0)
        self.W = KernelRidge(kernel='linear', alpha=1.0).fit(x, (U[:, :self.k]*S[:self.k]))
    def predict(self, x): return self.W.predict(x) @ self.Vk + self.mu
class _KNN:
    def __init__(self, k): self.k=k
    def fit(self, x, y): self.x=x; self.y=y
    def predict(self, x): return _knn_pred(self.x, self.y, x, self.k)
def _knn_pred(Xtr, Ytr, Xte, k):
    Xn = Xtr/ (np.linalg.norm(Xtr,axis=1,keepdims=True)+1e-9)
    Xt = Xte/ (np.linalg.norm(Xte,axis=1,keepdims=True)+1e-9)
    sim = Xt @ Xn.T
    out = np.zeros((len(Xte), Ytr.shape[1]), np.float32)
    for i in range(len(Xte)):
        nn = np.argsort(-sim[i])[:k]
        out[i] = Ytr[nn].mean(0)
    return out
results = {}
for dname, (X, Y) in directions.items():
    fold_ids = np.array([folds[t] for t in targets])
    oof = {}  # out-of-fold predictions per model
    for f in range(5):
        te = np.where(fold_ids==f)[0]; tr = np.where(fold_ids!=f)[0]
        preds = predict_models(X[tr], Y[tr], X[te])
        for mname, pr in preds.items():
            oof.setdefault(mname, {})[f] = (te, pr)
    res = {}
    # save OOF predictions for coverage calibration
    for mname, parts in oof.items():
        P = np.zeros(Y.shape, np.float32)
        for f, (te, pr) in parts.items(): P[te] = pr
        np.save(f'results/oof_{task}_' + dname.replace('/','_').replace('>','_').replace('+','_') + '__' + mname + '.npy', P)
    for mname, parts in oof.items():
        pearsons, r2s, top50, resid = [], [], [], []
        for f, (te, pr) in parts.items():
            for j, i in enumerate(te):
                o = Y[i]; p = pr[j]
                pearsons.append(float(np.corrcoef(p, o)[0,1]))
                cm = Y[[t2 for t2 in range(len(targets)) if folds[targets[t2]]!=f]].mean(0)
                r2s.append(float(1 - ((o-p)**2).sum()/max(((o-cm)**2).sum(),1e-9)))
                top50.append(len(set(np.argsort(-np.abs(o))[:50]) & set(np.argsort(-np.abs(p))[:50]))/50)
                resid.append(np.abs(o-p))
        res[mname] = {'pearson': pearsons, 'r2': r2s, 'top50': top50,
                      'resid_q90': float(np.quantile(np.concatenate(resid), 0.9))}
    results[dname] = res
    print(dname, 'done', time.time()-t0, 's', flush=True)
json.dump(results, open(f'results/r1_cv_{task}.json','w'))
print('TASK', task, 'DONE', time.time()-t0, 's', flush=True)
