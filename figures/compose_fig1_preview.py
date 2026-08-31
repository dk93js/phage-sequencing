"""
Layout preview ONLY (not the deliverable). Places 1A and the new _sized B/C
panels at 1:1 (native 300-DPI size) = the SAME zoom for all three. Because the
panels now share rcParams, their text is identical in pixels, so at equal zoom
the fonts are guaranteed to match. This is the exact condition to reproduce in
PowerPoint: place all three at the same zoom (e.g. 100%), do not stretch one
panel more than another. Deliverables are the standalone _sized PNGs.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

PANELS = Path(__file__).resolve().parent.parent / "panels"
DPI = 300

imgA = mpimg.imread(PANELS / "panelA_v3_smallmultiples.png")
imgB = mpimg.imread(PANELS / "cluster_trajectory_pooled_sized.png")
imgC = mpimg.imread(PANELS / "clonal_takeover_stack_pooled_sized.png")


def size_in(im):  # native (width, height) in inches at 300 DPI
    return im.shape[1] / DPI, im.shape[0] / DPI


aw, ah = size_in(imgA)
bw, bh = size_in(imgB)
cw, ch = size_in(imgC)

GAP_H = 0.45
GAP_V = 0.55
bottom_w = bw + GAP_H + cw
FIG_W = max(aw, bottom_w)
FIG_H = ah + GAP_V + max(bh, ch)

fig = plt.figure(figsize=(FIG_W, FIG_H))


def add(im, x0_in, y0_in, w_in, h_in):
    ax = fig.add_axes([x0_in / FIG_W, y0_in / FIG_H, w_in / FIG_W, h_in / FIG_H])
    ax.imshow(im)
    ax.axis("off")


add(imgA, (FIG_W - aw) / 2, FIG_H - ah, aw, ah)          # A across the top
bx = (FIG_W - bottom_w) / 2
add(imgB, bx, 0.0, bw, bh)                                # B | C below, 1:1
add(imgC, bx + bw + GAP_H, 0.0, cw, ch)

for letter, x_in, y_in in [("A", (FIG_W - aw) / 2 + 0.02, FIG_H - 0.02),
                           ("B", bx + 0.02, max(bh, ch) - 0.02),
                           ("C", bx + bw + GAP_H + 0.02, max(bh, ch) - 0.02)]:
    fig.text(x_in / FIG_W, y_in / FIG_H, letter, fontsize=20, fontweight="bold",
             va="top", ha="left", family="Arial")

fig.savefig(PANELS / "fig1_sized_preview.png", dpi=DPI, bbox_inches="tight")
print("wrote panels/fig1_sized_preview.png  (all panels at 1:1 / native 300 DPI)")
