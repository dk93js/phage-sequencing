"""R2 point 3 - the 215 cluster-discordant clonotypes: analyse them instead of
just excluding them from the per-cluster panels."""
import numpy as np, pandas as pd
from scipy import stats

cl = pd.read_csv("data/clusters_pooled.csv")
lab = cl.set_index(["seq", "lib"])["cluster"]
ppm = {r: cl.pivot_table(index="seq", columns="lib", values=f"PPM{r}") for r in range(5)}

shared = ppm[3].dropna().index                        # detected in both libraries
t = ppm[3].loc[shared]
d = pd.DataFrame(index=shared)
d["ss_lab"] = lab.reindex([(s, "ss") for s in shared]).values
d["ds_lab"] = lab.reindex([(s, "ds") for s in shared]).values
for r in range(5):
    p = ppm[r].loc[shared]
    d[f"lr{r}"] = np.log10(p["ss"] / p["ds"])
    d[f"ss{r}"], d[f"ds{r}"] = p["ss"], p["ds"]
d["disc"] = d.ss_lab != d.ds_lab
print(f"shared clonotypes: {len(d)}   discordant: {d.disc.sum()} ({d.disc.mean()*100:.1f}%)")

print("\nconfusion matrix (rows = dsDNA-copy cluster, cols = ssDNA-copy cluster):")
conf = pd.crosstab(d.ds_lab, d.ss_lab)
print(conf.to_string())

# ELISA calls
e = pd.read_excel("data/phage_elisa_wells.xlsx", sheet_name="ELISA")
e["Positivity"] = pd.to_numeric(e["Positivity"], errors="coerce")
e = e.dropna(subset=["VH CDR3", "Positivity"])
call = (e.groupby("VH CDR3")["Positivity"].mean() >= 0.5).astype(int)
d["elisa"] = call.reindex(d.index)

print("\nper transition class (R3 / R4 log10 ss/ds, median [IQR]):")
rows = []
for (ds_c, ss_c), g in d.groupby(["ds_lab", "ss_lab"]):
    if len(g) < 5: continue
    rows.append(dict(transition=f"ds C{ds_c} -> ss C{ss_c}", n=len(g),
                     kind="concordant" if ds_c == ss_c else "discordant",
                     lr_R3=round(g.lr3.median(), 3), lr_R4=round(g.lr4.median(), 3),
                     med_ssPPM_R3=round(g.ss3.median(), 1), med_dsPPM_R3=round(g.ds3.median(), 1),
                     n_elisa=int(g.elisa.notna().sum()), n_pos=int((g.elisa == 1).sum())))
tab = pd.DataFrame(rows).sort_values(["kind", "n"], ascending=[True, False])
print(tab.to_string(index=False))
tab.to_csv("out/r2_03_transition_classes.csv", index=False)

conc, disc = d[~d.disc], d[d.disc]
for r in (3, 4):
    U = stats.mannwhitneyu(disc[f"lr{r}"], conc[f"lr{r}"], alternative="two-sided")
    print(f"\nR{r}: median log10(ss/ds)  discordant {disc[f'lr{r}'].median():+.3f} "
          f"vs concordant {conc[f'lr{r}'].median():+.3f}   P = {U.pvalue:.2g}")

# does excluding them change the C1 result?
print("\nC1 slope at R3/R4 under three membership definitions:")
for name, mask in [("concordant C1 only (as published)", (d.ss_lab == 1) & (d.ds_lab == 1)),
                   ("any copy labelled C1",              (d.ss_lab == 1) | (d.ds_lab == 1)),
                   ("dsDNA-copy label only",              d.ds_lab == 1),
                   ("ssDNA-copy label only",              d.ss_lab == 1)]:
    out = []
    for r in (3, 4):
        g = d[mask & (d[f"ss{r}"] > 0) & (d[f"ds{r}"] > 0)]
        x, y = np.log10(g[f"ss{r}"]), np.log10(g[f"ds{r}"])
        out.append(f"R{r}: n={len(g)} slope={stats.linregress(x, y).slope:.3f} "
                   f"rho={stats.spearmanr(x, y).statistic:.3f}")
    print(f"  {name:36s} " + "   ".join(out))
d.to_csv("out/r2_03_shared_clonotype_table.csv")
