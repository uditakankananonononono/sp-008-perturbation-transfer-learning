import h5py, numpy as np, pandas as pd, json, time, resource, os
t0=time.time()
df = pd.read_pickle('/tmp/sp8_frangieh_obs.pkl')
targets = json.load(open('/tmp/sp8_frangieh_targets.json'))
conds = ['Control','IFNγ','Co-culture']
NG = 23712
base = (df.ngenes>=500)&(df.ncounts>=1000)&(df.mito<=20)
is_ctl = (df.pert=='control')
qc_targeting = base & (~is_ctl) & (df.nperts==1)
qc_ctl = base & is_ctl & (df.nperts==0)
tidx = {t:i for i,t in enumerate(targets)}; cidx = {c:i for i,c in enumerate(conds)}
pert = df['pert'].values; cond = df['cond'].values
row_T = np.full(len(df), -1, np.int32); row_C = np.full(len(df), -1, np.int32)
pt = pert; cd = cond
qct = qc_targeting.values; qcc = qc_ctl.values
for i in range(len(df)):
    if qcc[i] and cd[i] in cidx: row_C[i] = cidx[cd[i]]
    elif qct[i] and pt[i] in tidx and cd[i] in cidx: row_T[i] = tidx[pt[i]]*3 + cidx[cd[i]]
os.makedirs('results/accum', exist_ok=True)
T = np.memmap('results/accum/T_targetXcond.f32', dtype='float32', mode='w+', shape=(len(targets)*3, NG))
C = np.zeros((3, NG), np.float32)
Tn = np.bincount(row_T[row_T>=0], minlength=len(targets)*3).astype(np.int32)
Cn = np.bincount(row_C[row_C>=0], minlength=3).astype(np.int32)
with h5py.File('data/raw/FrangiehIzar2021_RNA.h5ad','r') as f:
    X = f['X']; data = X['data']; indices = X['indices']; indptr = X['indptr']
    for j in range(NG):
        i0, i1 = int(indptr[j]), int(indptr[j+1])
        if i0 == i1: continue
        rows = indices[i0:i1]; vals = data[i0:i1]
        gT = row_T[rows]; m = gT>=0
        if m.any(): np.add.at(T[:, j], gT[m], vals[m])
        gC = row_C[rows]; m2 = gC>=0
        if m2.any(): np.add.at(C[:, j], gC[m2], vals[m2])
        if j % 4000 == 0: print(f'col {j}/{NG} {time.time()-t0:.0f}s RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f}GB', flush=True)
T.flush()
np.save('results/accum/Tn.npy', Tn); np.save('results/accum/C.npy', C); np.save('results/accum/Cn.npy', Cn)
json.dump({'targets': targets, 'conditions': conds}, open('results/accum/index.json','w'))
print('DONE', time.time()-t0, 's; peak RSS GB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6, flush=True)
