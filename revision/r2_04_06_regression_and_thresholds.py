"""R2 point 4 - symmetric (errors-in-both-variables) regression instead of OLS.
R2 point 6 - sensitivity of the conclusions to the PPM abundance cut-off."""
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster

cl = pd.read_csv("data/clusters_pooled.csv")

def sma(x, y):
    r = np.corrcoef(x, y)[0, 1]
    return np.sign(r) * y.std(ddof=1) / x.std(ddof=1), r

def deming(x, y):                                    # orthogonal regression, delta = 1
    sxx, syy = x.var(ddof=1), y.var(ddof=1)
    sxy = np.cov(x, y, ddof=1)[0, 1]
    return ((syy - sxx) + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)

# ---------- R2-4: three estimators, published cluster labels -----------------
lab = cl.set_index(["seq", "lib"])["cluster"]
ppm = {r: cl.pivot_table(index="seq", columns="lib", values=f"PPM{r}") for r in range(5)}
rows = []
for r in range(5):
    t = ppm[r].dropna(); t = t[(t.ss > 0) & (t.ds > 0)]
    ls = lab.reindex([(s, "ss") for s in t.index]).values
    ld = lab.reindex([(s, "ds") for s in t.index]).values
    for name, m in [("all", np.ones(len(t), bool)), ("C1", (ls == ld) & (ls == 1)),
                    ("C2", (ls == ld) & (ls == 2)), ("C3", (ls == ld) & (ls == 3))]:
        x, y = np.log10(t.ss.values[m]), np.log10(t.ds.values[m])
        s_sma, pear = sma(x, y)
        rows.append(dict(cluster=name, round=f"R{r}", n=int(m.sum()),
                         spearman_rho=round(stats.spearmanr(x, y).statistic, 3),
                         pearson_r=round(pear, 3),
                         OLS=round(stats.linregress(x, y).slope, 3),
                         SMA=round(s_sma, 3), Deming=round(deming(x, y), 3)))
reg = pd.DataFrame(rows)
print("R2-4  slope estimators (OLS as published vs symmetric alternatives)\n")
print(reg.pivot(index="cluster", columns="round",
                values=["OLS", "SMA", "Deming"]).to_string())
reg.to_csv("out/r2_04_slope_estimators.csv", index=False)

# ---------- R2-6: abundance-threshold sensitivity ----------------------------
print("\n\nR2-6  PPM cut-off sensitivity (filter -> re-cluster -> re-measure)\n")
P = cl[[f"PPM{r}" for r in range(5)]]
out = []
for thr in (100, 200, 500, 1000):
    d = cl[P.max(axis=1) >= thr].reset_index(drop=True)
    M = d[[f"PPM{r}" for r in range(5)]].values
    Zs = (M - M.mean(1, keepdims=True)) / M.std(1, keepdims=True)
    k3 = fcluster(linkage(Zs, method="ward"), 3, criterion="maxclust")
    gain = {c: (Zs[k3 == c][:, 4] - Zs[k3 == c][:, 0]).mean() for c in np.unique(k3)}
    enr = max(gain, key=gain.get)
    d["enr"] = k3 == enr
    lm = d.set_index(["seq", "lib"])["enr"]
    pv = {r: d.pivot_table(index="seq", columns="lib", values=f"PPM{r}") for r in range(5)}
    for r in (2, 3, 4):
        t = pv[r].dropna(); t = t[(t.ss > 0) & (t.ds > 0)]
        m = (lm.reindex([(s, "ss") for s in t.index]).values &
             lm.reindex([(s, "ds") for s in t.index]).values)
        x, y = np.log10(t.ss.values[m]), np.log10(t.ds.values[m])
        xa, ya = np.log10(t.ss.values), np.log10(t.ds.values)
        out.append(dict(PPM_cut=thr, rows=len(d), round=f"R{r}", n_enriching=int(m.sum()),
                        slope_enriching=round(stats.linregress(x, y).slope, 3),
                        rho_enriching=round(stats.spearmanr(x, y).statistic, 3),
                        slope_all=round(stats.linregress(xa, ya).slope, 3),
                        rho_all=round(stats.spearmanr(xa, ya).statistic, 3)))
th = pd.DataFrame(out)
print(th.to_string(index=False))
th.to_csv("out/r2_06_threshold_sensitivity.csv", index=False)
