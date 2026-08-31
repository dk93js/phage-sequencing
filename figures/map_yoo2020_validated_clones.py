"""
Map the experimentally-validated AR binders from Yoo et al. 2020 Table 2 onto our
HiSeq CDRH3 ss/ds analysis set: are the validated binders captured, what cluster do
they fall in, and are they themselves subject to the ssDNA/dsDNA bias?

Yoo 2020 Table 2 unique HCDR3 AA sequences (AR1-AR15; AR5-AR15 all share one HCDR3):
    AR1               GSGGVDSIDA
    AR2               SADGYGWDTAGNMDA
    AR3               TAGTCTTSCNAGAYIDA
    AR4               TTCSGSYGWCADSIDA
    AR5-AR15 (11/15)  SADSCATCATYPSEIDT   <- dominant validated binder
Also probe SADSCATCATYPSEIDA, the high-abundance 1-residue (A vs T) variant currently
used as the "unbiased" worked clone.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

CODE_DIR = Path(__file__).resolve().parent
DATA = CODE_DIR.parent / "data"
PANELS = CODE_DIR.parent / "panels"
CL = {1: "C1 enriching", 2: "C2 intermediate", 3: "C3 depleted"}
C_CLUSTER = {1: "#E63946", 2: "#457B9D", 3: "#2A9D8F"}

mpl.rcParams.update({
    "font.family": "Arial", "font.size": 13, "axes.titlesize": 15,
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

VALIDATED = {
    "GSGGVDSIDA": "AR1",
    "SADGYGWDTAGNMDA": "AR2",
    "TAGTCTTSCNAGAYIDA": "AR3",
    "TTCSGSYGWCADSIDA": "AR4",
    "SADSCATCATYPSEIDT": "AR5-15 (dominant)",
}
VARIANT = {"SADSCATCATYPSEIDA": "1-aa variant of AR5-15 (A vs T); not in Table 2"}

p = pd.read_csv(DATA / "clusters_pooled.csv")

def describe(seq, tag):
    sub = p[p.seq == seq]
    if len(sub) == 0:
        print(f"  {seq:20s} [{tag}]  -> NOT in filtered analysis set (q40, PPM>=100)")
        return None
    rows = {r.lib: r for _, r in sub.iterrows()}
    ss_c = int(rows["ss"].cluster) if "ss" in rows else None
    ds_c = int(rows["ds"].cluster) if "ds" in rows else None
    ss3 = rows["ss"].PPM3 if "ss" in rows else 0
    ds3 = rows["ds"].PPM3 if "ds" in rows else 0
    libs = "+".join(rows.keys())
    fold = (f"{max(ss3,ds3)/max(min(ss3,ds3),1):.1f}x "
            f"{'ss' if ss3>ds3 else 'ds'}-leaning") if ss3 and ds3 else "single-lib"
    print(f"  {seq:20s} [{tag}]")
    print(f"      libs={libs}  cluster ss={CL.get(ss_c,'-')} / ds={CL.get(ds_c,'-')}"
          f"  | R3 PPM ss={int(ss3)} ds={int(ds3)}  ({fold})")
    return {"seq": seq, "tag": tag, "libs": libs, "ss_cluster": ss_c,
            "ds_cluster": ds_c, "ss_PPM3": int(ss3), "ds_PPM3": int(ds3)}

print("=== Yoo 2020 Table 2 validated binders in our analysis set ===")
recs = [r for seq, tag in VALIDATED.items() if (r := describe(seq, tag))]
print("\n=== high-abundance 1-residue variant ===")
recs += [r for seq, tag in VARIANT.items() if (r := describe(seq, tag))]

found = sum(1 for seq in VALIDATED if len(p[p.seq == seq]))
print(f"\nCaptured: {found}/{len(VALIDATED)} unique Table-2 HCDR3s pass our filter.")
if recs:
    out = pd.DataFrame(recs)
    out.to_csv(DATA / "yoo2020_table2_clone_mapping.csv", index=False)
    print(f"wrote {DATA / 'yoo2020_table2_clone_mapping.csv'}")


# =========================================================
# Panel: validated binders overlaid on the full ss-vs-ds map (R3)
# =========================================================
def render_validated_overlay(round_=3):
    piv = p.pivot_table(index="seq", columns="lib", values=f"PPM{round_}",
                        aggfunc="first").fillna(0)
    clu = p.pivot_table(index="seq", columns="lib", values="cluster",
                        aggfunc="first")
    both = piv[(piv["ss"] > 0) & (piv["ds"] > 0)]

    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(both["ss"], both["ds"], s=10, color="#cfcfcf", alpha=0.6,
               edgecolor="none", zorder=1, label=f"all shared clonotypes (n={len(both)})")
    lo, hi = 50, max(both.max().max(), 2e5) * 1.5
    ax.plot([lo, hi], [lo, hi], color="#555", linewidth=1.1, ls="-", zorder=2,
            label="y = x")

    # stars for the 5 validated binders, colored by (concordant) cluster
    label_offsets = {  # hand-nudged so labels don't collide
        "AR1": (1.25, 0.78), "AR2": (1.3, 1.0), "AR3": (0.42, 1.15),
        "AR4": (1.3, 0.72), "AR5-15 (dominant)": (0.30, 1.35),
    }
    for seq, tag in VALIDATED.items():
        if seq not in piv.index:
            continue
        x, y = piv.loc[seq, "ss"], piv.loc[seq, "ds"]
        c = int(clu.loc[seq, "ss"])
        ax.scatter([x], [y], s=240, marker="*", color=C_CLUSTER[c],
                   edgecolor="black", linewidth=0.8, zorder=4)
        fx, fy = label_offsets.get(tag, (1.2, 1.2))
        ax.annotate(tag.replace(" (dominant)", "*"), (x, y), (x * fx, y * fy),
                    fontsize=9, fontweight="bold", zorder=5,
                    arrowprops=dict(arrowstyle="-", color="#888", lw=0.6))
    # the high-abundance 1-aa variant of the dominant binder
    if "SADSCATCATYPSEIDA" in piv.index:
        x, y = piv.loc["SADSCATCATYPSEIDA", "ss"], piv.loc["SADSCATCATYPSEIDA", "ds"]
        ax.scatter([x], [y], s=240, marker="*", color=C_CLUSTER[1],
                   edgecolor="black", linewidth=0.8, zorder=4)
        ax.annotate("AR5-15 -A variant\n(deep-seq dominant)", (x, y),
                    (x * 0.12, y * 0.55), fontsize=9, fontweight="bold", zorder=5,
                    arrowprops=dict(arrowstyle="-", color="#888", lw=0.6))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("ssDNA PPM")
    ax.set_ylabel("dsDNA PPM")
    ax.set_title(f"Yoo 2020 validated binders on the ss-vs-ds map (R{round_})\n"
                 "stars colored by cluster; * = AR5-15 dominant HCDR3", fontsize=12)
    handles = [plt.Line2D([], [], marker="*", color=C_CLUSTER[k], ls="",
                          markeredgecolor="black", markersize=13, label=CL[k])
               for k in (1, 2, 3)]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=10)
    fig.savefig(PANELS / "yoo2020_validated_binders_scatter_R3.png")
    plt.close(fig)
    print(f"wrote {PANELS / 'yoo2020_validated_binders_scatter_R3.png'}")


render_validated_overlay(3)
