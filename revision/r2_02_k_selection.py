"""R2 point 2 - why k = 3 for the Ward clustering?

(a) internal validity indices across k = 2..8
(b) gap statistic (Tibshirani, PCA-aligned uniform reference)
(c) bootstrap cluster stability (Hennig's clusterboot Jaccard)
(d) does the paper's conclusion survive a different k?
"""
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy import stats
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score

rng = np.random.default_rng(0)
cl = pd.read_csv("data/clusters_pooled.csv")
Z = cl[[f"z{r}" for r in range(5)]].values          # z-scored 5-round trajectories
P = cl[[f"PPM{r}" for r in range(5)]].values
n = len(Z)
print(f"n = {n} clonotype-by-library rows, {Z.shape[1]} rounds")

def ward(X, k):
    return fcluster(linkage(X, method="ward"), k, criterion="maxclust")

# ---- (a) internal indices ---------------------------------------------------
D = pdist(Z)
rows = []
for k in range(2, 9):
    lab = ward(Z, k)
    rows.append(dict(k=k,
                     silhouette=round(silhouette_score(Z, lab), 4),
                     calinski_harabasz=round(calinski_harabasz_score(Z, lab), 1),
                     davies_bouldin=round(davies_bouldin_score(Z, lab), 4),
                     min_cluster_size=int(np.bincount(lab)[1:].min())))
idx = pd.DataFrame(rows)

# ---- (b) gap statistic ------------------------------------------------------
def Wk(X, lab):
    return sum(((X[lab == c] - X[lab == c].mean(0))**2).sum() for c in np.unique(lab))

Xc = Z - Z.mean(0)
U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
Xp = Xc @ Vt.T
lo, hi = Xp.min(0), Xp.max(0)
B = 50
gaps, sks = [], []
for k in range(2, 9):
    logW = np.log(Wk(Z, ward(Z, k)))
    ref = []
    for _ in range(B):
        Rp = rng.uniform(lo, hi, size=Xp.shape)
        R = Rp @ Vt + Z.mean(0)
        ref.append(np.log(Wk(R, ward(R, k))))
    ref = np.array(ref)
    gaps.append(ref.mean() - logW)
    sks.append(ref.std() * np.sqrt(1 + 1 / B))
idx["gap"] = np.round(gaps, 4)
idx["gap_se"] = np.round(sks, 4)
# Tibshirani rule: smallest k with gap(k) >= gap(k+1) - s(k+1)
pick = next((k for i, k in enumerate(range(2, 9))
             if i + 1 < len(gaps) and gaps[i] >= gaps[i + 1] - sks[i + 1]), None)
print(f"\ngap-statistic choice (Tibshirani 1SE rule): k = {pick}")
print(idx.to_string(index=False))
idx.to_csv("out/r2_02_k_indices.csv", index=False)

# ---- (c) bootstrap stability ------------------------------------------------
B2 = 200
jac = {k: np.zeros((B2, k)) for k in range(2, 7)}
for b in range(B2):
    samp = rng.integers(0, n, n)
    Zb = Z[samp]
    for k in range(2, 7):
        base = ward(Z, k)[samp]                    # original labels of the resampled points
        new = ward(Zb, k)
        for c in range(1, k + 1):
            A = base == c
            best = max(((A & (new == d)).sum() / max((A | (new == d)).sum(), 1))
                       for d in range(1, k + 1))
            jac[k][b, c - 1] = best
print("\nbootstrap cluster stability (mean Jaccard over 200 resamples; >0.75 = stable):")
stab = []
for k in range(2, 7):
    m = jac[k].mean(0)
    stab.append(dict(k=k, mean_jaccard=round(m.mean(), 3),
                     per_cluster=" / ".join(f"{v:.3f}" for v in np.sort(m)[::-1]),
                     min_cluster=round(m.min(), 3)))
stabdf = pd.DataFrame(stab)
print(stabdf.to_string(index=False))
stabdf.to_csv("out/r2_02_bootstrap_stability.csv", index=False)

# ---- (d) does the conclusion survive another k? -----------------------------
ppm = {r: cl.pivot_table(index="seq", columns="lib", values=f"PPM{r}") for r in range(5)}
print("\nkey result (slope of the most-enriching cluster) as a function of k:")
res = []
for k in range(2, 7):
    lab = ward(Z, k)
    d = cl[["seq", "lib"]].copy(); d["cl"] = lab
    gain = pd.Series({c: (Z[lab == c][:, 4] - Z[lab == c][:, 0]).mean() for c in np.unique(lab)})
    top = gain.idxmax()
    d["enr"] = d["cl"] == top
    lm = d.set_index(["seq", "lib"])["enr"]
    for r in (2, 3, 4):
        t = ppm[r].dropna(); t = t[(t["ss"] > 0) & (t["ds"] > 0)]
        ss_l = lm.reindex([(s, "ss") for s in t.index]).values
        ds_l = lm.reindex([(s, "ds") for s in t.index]).values
        m = ss_l & ds_l
        x, y = np.log10(t["ss"].values[m]), np.log10(t["ds"].values[m])
        res.append(dict(k=k, round=f"R{r}", n=int(m.sum()),
                        slope=round(stats.linregress(x, y).slope, 3),
                        rho=round(stats.spearmanr(x, y).statistic, 3)))
resdf = pd.DataFrame(res).pivot(index="k", columns="round", values=["n", "slope", "rho"])
print(resdf.to_string())
pd.DataFrame(res).to_csv("out/r2_02_conclusion_vs_k.csv", index=False)
