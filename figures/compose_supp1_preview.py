"""
Layout preview ONLY (not the deliverable) for Supplementary Figure 1 panels
B (3 enrichment heatmaps) and C (2 per-library cluster trajectories, _sized).
Panel A (the manual stats table) is not included. All panels at 1:1 (native
300 DPI) so apparent font sizes are faithful.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PANELS = Path(__file__).resolve().parent.parent / "panels"
DPI = 300

hm = [mpimg.imread(PANELS / f"enrichment_heatmap_{d}.png") for d in ["ssDNA", "dsDNA", "pooled"]]
tr = [mpimg.imread(PANELS / f"cluster_trajectory_{d}_sized.png") for d in ["ssDNA", "dsDNA"]]


def size_in(im):
    return im.shape[1] / DPI, im.shape[0] / DPI


GAP_H = 0.30
GAP_V = 0.55
hm_sz = [size_in(im) for im in hm]
tr_sz = [size_in(im) for im in tr]
b_w = sum(w for w, h in hm_sz) + GAP_H * (len(hm) - 1)
b_h = max(h for w, h in hm_sz)
c_w = sum(w for w, h in tr_sz) + GAP_H * (len(tr) - 1)
c_h = max(h for w, h in tr_sz)

FIG_W = max(b_w, c_w)
FIG_H = b_h + GAP_V + c_h
fig = plt.figure(figsize=(FIG_W, FIG_H))


def add(im, x0, y0, w, h):
    ax = fig.add_axes([x0 / FIG_W, y0 / FIG_H, w / FIG_W, h / FIG_H])
    ax.imshow(im)
    ax.axis("off")


bx = (FIG_W - b_w) / 2
for im, (w, h) in zip(hm, hm_sz):
    add(im, bx, FIG_H - h, w, h)          # B row, top-aligned
    bx += w + GAP_H
cx = (FIG_W - c_w) / 2
for im, (w, h) in zip(tr, tr_sz):
    add(im, cx, 0.0, w, h)                # C row, centered below
    cx += w + GAP_H

fig.text((FIG_W - b_w) / 2 / FIG_W, (FIG_H - 0.02) / FIG_H, "B",
         fontsize=22, fontweight="bold", va="top", ha="left", family="Arial")
fig.text((FIG_W - c_w) / 2 / FIG_W, (c_h - 0.02) / FIG_H, "C",
         fontsize=22, fontweight="bold", va="top", ha="left", family="Arial")

fig.savefig(PANELS / "supp_fig1_BC_preview.png", dpi=200, bbox_inches="tight")
print("wrote panels/supp_fig1_BC_preview.png")
print(f"trajectory dims: ss={tr[0].shape[1]}x{tr[0].shape[0]}px  ds={tr[1].shape[1]}x{tr[1].shape[0]}px")
print(f"heatmap heights: {[im.shape[0] for im in hm]}px")
