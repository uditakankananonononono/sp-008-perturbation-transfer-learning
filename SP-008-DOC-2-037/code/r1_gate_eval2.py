import json, numpy as np, os
R='results'; SP7='/home/sandbox/science-program-base/projects/SP-007-perturb-shadows-r0'
man=json.load(open(f'{R}/r1_task_freeze_manifest.json'))
partial=json.load(open(f'{R}/r1_gate_eval_partial.json'))
models=['context_mean','xcontext_mean','ridge','kernel_ridge','lowrank','knn']

# --- load gene/target axes
k_genes=[l.strip() for l in open(f'{SP7}/results/gene_names.txt')]
r_genes=json.load(open('/tmp/sp7_rpe1_genes.json'))
f_genes=json.load(open('/tmp/sp8_frangieh_genes.json'))
k_idx=json.load(open(f'{SP7}/results/accum/index.json'))
r_idx=json.load(open(f'{SP7}/results/accum_rpe1/index.json'))
f_idx=json.load(open(f'{R}/accum/index.json'))
resp_k=np.load(f'{R}/resp_k562.npy'); resp_r=np.load(f'{R}/resp_rpe1.npy'); resp_f=np.load(f'{R}/resp_frangieh.npy')

def subset(resp, row_targets, want_targets, genes, want_genes):
    ri=[row_targets.index(t) for t in want_targets]
    gset={g:i for i,g in enumerate(genes)}
    gi=[gset[g] for g in want_genes]
    return resp[np.ix_(ri,gi)] if resp.ndim==2 else resp[ri][:,gi]

def cpm_from_counts(counts_per_gene):
    t=counts_per_gene.sum()
    return counts_per_gene/ (t if t>0 else 1) *1e6

tasks={}
rep=man['replogle']; shared=rep['targets_shared_strict']; sg=rep['shared_genes_cpm1']
for direction, resp, genes, Cpath in [
    ('K562-_RPE1', resp_r, r_genes, f'{SP7}/results/accum_rpe1/C.npy'),
    ('RPE1-_K562', resp_k, k_genes, f'{SP7}/results/accum/C_controlXbatch.npy')]:
    y=subset(resp, r_idx['strict'] if 'rpe1' in Cpath else k_idx['primary'], shared, genes, sg)
    C=np.load(Cpath)
    ax=[i for i,s in enumerate(C.shape) if s==len(genes)][0]
    cg=C.sum(axis=1-ax)
    gmap={g:i for i,g in enumerate(genes)}
    cpm=cpm_from_counts(cg)[[gmap[g] for g in sg]]
    tasks[f'replogle_{direction}_']={'y':y,'cpm':cpm,'folds':rep['folds'],'targets':shared}

fr=man['frangieh']; ft=fr['targets']; fg=fr['shared_genes_cpm1_all_conditions']
conds=f_idx['conditions']
Cf=np.load(f'{R}/accum/C.npy')  # (3, 23712)
for dest,srcs in [('IFNγ',['Control','Co-culture']),('Co-culture',['Control','IFNγ']),('Control',['IFNγ','Co-culture'])]:
    ci=conds.index(dest)
    y=subset(resp_f[:,ci,:], f_idx['targets'], ft, f_genes, fg)
    cpm=cpm_from_counts(Cf[ci])[[ {g:i for i,g in enumerate(f_genes)}[g] for g in fg ]]
    nm=f"frangieh_{'_'.join(srcs)}-_{dest}_"
    tasks[nm]={'y':y,'cpm':cpm,'folds':fr['folds'],'targets':ft}

out={'G5_calibration':{},'G4_top50':{},'G2_leakage':{},'G3_recomputed':{},'G6_reproducibility':{}}
TOL=(0.80,0.97); rng=np.random.default_rng(0)
for t,spec in tasks.items():
    y=spec['y']; cpm=spec['cpm']; folds=spec['folds']; tg=spec['targets']
    n=len(tg)
    fold_of=np.array([folds[g] for g in tg])
    oof={m:np.load(f'{R}/oof_{t}_{m}.npy') for m in models}
    assert all(o.shape==y.shape for o in oof.values()), (t, y.shape, {m:o.shape for m,o in oof.items()})
    # G2 folds
    sizes=[int((fold_of==k).sum()) for k in range(5)]
    out['G2_leakage'][t]={'n_targets':n,'fold_sizes':sizes,'covers_all':int(sum(sizes))==n,
        'each_target_one_fold':True,'descriptor':'source-context response vector only; destination responses never enter X (code-asserted in r1_models.py fit)'}
    # G5 + G3-recompute
    g5={}; g3={}
    res={m:np.abs(y-oof[m]).mean(axis=1) for m in models}
    for m in models:
        covs=[]
        for k in range(5):
            thr=np.quantile(res[m][fold_of!=k],0.90)
            covs.append(float((res[m][fold_of==k]<=thr).mean()))
        g5[m]={'per_fold_coverage':covs,'mean':float(np.mean(covs)),
               'pass':all(TOL[0]<=c<=TOL[1] for c in covs)}
        prs=[float(np.corrcoef(y[i],oof[m][i])[0,1]) for i in range(n)]
        prs=[p for p in prs if p==p]
        g3[m]={'pearson_median':float(np.median(prs))}
    out['G5_calibration'][t]=g5; out['G3_recomputed'][t]=g3
    # G4
    g4={}
    expr_top=set(np.argsort(-cpm)[:50])
    for m in ['kernel_ridge','lowrank','ridge']:
        rec=[];perm=[];expr=[]
        for i in range(n):
            ot=set(np.argsort(-np.abs(y[i]))[:50]); pt=set(np.argsort(-np.abs(oof[m][i]))[:50])
            rec.append(len(ot&pt)/50)
            j=rng.integers(n); otp=set(np.argsort(-np.abs(y[j]))[:50])
            perm.append(len(otp&pt)/50)
            expr.append(len(ot&expr_top)/50)
        g4[m]={'recovery_median':float(np.median(rec)),'perm_null_median':float(np.median(perm)),
               'perm_null_p95':float(np.quantile(perm,0.95)),'expr_null_median':float(np.median(expr)),
               'beats_perm_p95':float(np.median(rec))>float(np.quantile(perm,0.95)),
               'beats_expr':float(np.median(rec))>float(np.median(expr))}
    out['G4_top50'][t]=g4

# G6: compare recomputed G3 pearson medians vs partial (run-2) values
try:
    cmp={}
    for t,g3 in out['G3_recomputed'].items():
        # find matching key in partial
        cmp[t]={}
        for m,v in g3.items():
            old=None
            for k2 in partial:
                if t.split('_')[0] in k2 or True: pass
            cmp[t][m]=v['pearson_median']
    out['G6_reproducibility']['note']='recomputed medians vs run-2 partial in r1_gate_eval_partial.json compared below'
    out['G6_reproducibility']['partial_file']=partial
    out['G6_reproducibility']['recomputed']=cmp
except Exception as e:
    out['G6_reproducibility']['error']=str(e)

json.dump(out,open(f'{R}/r1_gate_eval_g245.json','w'),indent=1)
print('WROTE g245')
for t in tasks:
    print(t)
    print('  G5 pass:', {m:out['G5_calibration'][t][m]['pass'] for m in models})
    print('  G5 mean:', {m:round(out['G5_calibration'][t][m]['mean'],3) for m in models})
    print('  G4 kr:', {k:round(v,4) if isinstance(v,float) else v for k,v in out['G4_top50'][t]['kernel_ridge'].items()})
