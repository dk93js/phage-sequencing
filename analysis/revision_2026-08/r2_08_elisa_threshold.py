"""R2 point 8 - ELISA call rule: correct the '>=0.7' text error and show the
antigen-reactivity result is unchanged under every reasonable calling rule.

Reproduces the paper's numbers (mean call >= 0.5 -> 90 reactive / 478 non-reactive;
146 matched -> 35 +/111 -) and re-runs Figure 3B under alternative rules.
"""
import os
os.makedirs("out/figures", exist_ok=True)
import numpy as np, pandas as pd
from scipy import stats

ELISA = "data/phage_elisa_wells.xlsx"
CL    = "data/clusters_pooled.csv"

e = pd.read_excel(ELISA, sheet_name="ELISA")
e = e.dropna(subset=["VH CDR3", "Positivity"])
e["Positivity"] = pd.to_numeric(e["Positivity"], errors="coerce")
e = e.dropna(subset=["Positivity"])
print(f"assayed wells (rows with a call): {len(e)}   unique HCDR3: {e['VH CDR3'].nunique()}")

g = e.groupby("VH CDR3")
summ = pd.DataFrame({
    "n_wells":  g.size(),
    "mean":     g["Positivity"].mean(),
    "wfrac":    g.apply(lambda d: np.average(d["Positivity"], weights=d["Frequency"])),
    "dominant": g.apply(lambda d: d.loc[d["Frequency"].idxmax(), "Positivity"]),
})

rules = {
    "mean >= 0.5 (simple majority; used in all analyses)": summ["mean"] >= 0.5,
    "mean >= 0.7 (as mis-stated in the submitted Methods)": summ["mean"] >= 0.7,
    "frequency-weighted majority (>= 0.5)":                 summ["wfrac"] >= 0.5,
    "dominant-well call":                                   summ["dominant"] >= 0.5,
    "ties (mean == 0.5) excluded":                          summ["mean"] > 0.5,
}
drop_ties = summ["mean"] == 0.5

# --- paired ss/ds abundances -------------------------------------------------
cl = pd.read_csv(CL)
ppm = {r: cl.pivot_table(index="seq", columns="lib", values=f"PPM{r}") for r in range(5)}

rows = []
for name, call in rules.items():
    lab = call.astype(int)
    if name.startswith("ties"):
        lab = lab[~drop_ties]
    for r in (3, 4):
        t = ppm[r].dropna()
        t = t[(t["ss"] > 0) & (t["ds"] > 0)]
        idx = t.index.intersection(lab.index)
        lr = np.log10(t.loc[idx, "ss"] / t.loc[idx, "ds"])
        pos, neg = lr[lab.loc[idx] == 1], lr[lab.loc[idx] == 0]
        U = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        rows.append(dict(rule=name, round=f"R{r}", n_reactive_total=int(call.sum()),
                         n_matched=len(idx), n_pos=len(pos), n_neg=len(neg),
                         median_pos=round(pos.median(), 3), median_neg=round(neg.median(), 3),
                         P=f"{U.pvalue:.2g}"))

out = pd.DataFrame(rows)
print()
print(out.to_string(index=False))
out.to_csv("out/r2_08_elisa_rule_sensitivity.csv", index=False)

n_ties = int(drop_ties.sum())
print(f"\nreplicated HCDR3 (>=2 wells): {(summ['n_wells'] >= 2).sum()}   singletons: {(summ['n_wells'] == 1).sum()}")
print(f"exact 50:50 ties: {n_ties}  (called reactive by the >=0.5 rule)")
print("\ndisagreement between the >=0.5 and >=0.7 rules:")
diff = summ[(summ['mean'] >= 0.5) != (summ['mean'] >= 0.7)]
print(diff.assign(seq=diff.index).to_string(index=False))
