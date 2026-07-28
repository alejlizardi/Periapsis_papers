#!/usr/bin/env python3
"""Generate the gadget figures (E12) in the style of the union_viz.html teaching tool
(misc/PC_Ham_Teaching/tools/union_viz.html in the Periapsis repo).

Conventions (matching the tool and the paper's terminology):
  - vertices = positions 0..b-1 as dark dots on a gray baseline;
  - BLUE arcs above the line = P-edges {i,i+1};
  - RED arcs below the line = V-edges, arc width proportional to the position jump;
  - PURPLE paired arcs (one above, one below) = doubled pairs;
  - general-b panels elide the middle of the doubled run with a centered ellipsis,
    dashed purple stubs indicating the run continues through it.

Outputs: fig_gadget_Tb.pdf  (T_8 concrete + general even b)
         fig_gadget_Tpb.pdf (T'_9 concrete + general odd b)
Run from this directory:  python gen_gadget_figs.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

BLUE, RED, DOUBLE, INK, MUTED, LINE = "#2196f3", "#e53935", "#8e24aa", "#1a1a1a", "#666666", "#333333"
VSCALE = 0.30  # arch peak = VSCALE * span (tool default 0.42; reduced for print)


def arch(ax, xa, xb, y0, direction, color, lw, dashed=False, vscale=VSCALE):
    """Cubic-Bezier arch from (xa,y0) to (xb,y0); direction +1 up, -1 down."""
    h = direction * (4.0 / 3.0) * vscale * abs(xb - xa)
    path = Path([(xa, y0), (xa, y0 + h), (xb, y0 + h), (xb, y0)],
                [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, fill=False, edgecolor=color, lw=lw,
                           capstyle="round", alpha=0.95,
                           linestyle=(0, (2.4, 2.4)) if dashed else "solid"))


def panel(ax, slots, p_edges, v_edges, doubled, dashed_doubled=(), ellipsis_slot=None,
          value_labels=None, vscale=VSCALE):
    """Draw one H_sigma panel.

    slots: list of position-label strings ('' allowed), one per x-slot.
    p_edges / v_edges: (i, j) slot pairs for P-only (top, blue) / V-only (bottom, red).
    doubled: (i, j) slot pairs drawn as a purple arc pair (top + bottom).
    dashed_doubled: like doubled but dashed (continuation into the ellipsis).
    ellipsis_slot: slot index that holds a centered ellipsis instead of a vertex.
    value_labels: optional list of value-label strings under each vertex.
    """
    n = len(slots)
    top = max((abs(b - a) for a, b in list(p_edges) + list(doubled) + list(dashed_doubled)), default=1)
    bot = max((abs(b - a) for a, b in list(v_edges) + list(doubled) + list(dashed_doubled)), default=1)
    ax.set_xlim(-0.7, n - 0.3)
    ax.set_ylim(-vscale * bot - 0.55, vscale * top + 0.55)
    ax.set_aspect("equal")
    ax.axis("off")

    ax.plot([-0.35, n - 0.65], [0, 0], color="#bbbbbb", lw=0.8, zorder=1)

    for a, b in p_edges:
        arch(ax, a, b, 0, +1, BLUE, 1.6)
    for a, b in v_edges:
        arch(ax, a, b, 0, -1, RED, 1.5)
    for a, b in doubled:
        arch(ax, a, b, 0, +1, DOUBLE, 2.0)
        arch(ax, a, b, 0, -1, DOUBLE, 2.0)
    for a, b in dashed_doubled:
        arch(ax, a, b, 0, +1, DOUBLE, 1.6, dashed=True)
        arch(ax, a, b, 0, -1, DOUBLE, 1.6, dashed=True)

    for i, lab in enumerate(slots):
        if i == ellipsis_slot:
            ax.text(i, 0, r"$\cdots$", ha="center", va="center", fontsize=11, color=INK, zorder=3)
            continue
        ax.plot([i], [0], "o", ms=3.6, color=LINE, zorder=3)
        ax.text(i, vscale * top + 0.34, lab, ha="center", va="bottom",
                fontsize=9, color=INK, zorder=3)
        if value_labels and value_labels[i]:
            ax.text(i, -vscale * bot - 0.30, value_labels[i], ha="center", va="top",
                    fontsize=7, color=MUTED, zorder=3)


def concrete_edges(sigma):
    """P-only, V-only, doubled slot pairs for a concrete permutation (tool logic)."""
    n = len(sigma)
    inv = [0] * n
    for i, v in enumerate(sigma):
        inv[v] = i
    dbl = {(i, i + 1) for i in range(n - 1) if abs(sigma[i + 1] - sigma[i]) == 1}
    p_only = [(i, i + 1) for i in range(n - 1) if (i, i + 1) not in dbl]
    v_only = []
    for k in range(n - 1):
        a, b = sorted((inv[k], inv[k + 1]))
        if (a, b) not in dbl:
            v_only.append((a, b))
    return p_only, sorted(set(v_only)), sorted(dbl)


def fig_Tb():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 3.6))

    # (a) concrete T_8 = (0,7,2,3,4,5,6,1)
    T8 = [0, 7, 2, 3, 4, 5, 6, 1]
    p, v, d = concrete_edges(T8)
    panel(ax1, [str(i) for i in range(8)], p, v, d,
          value_labels=[str(x) for x in T8])

    # (b) general even b: slots 0..4, ellipsis, b-4..b-1
    labs = ["0", "1", "2", "3", "4", "", r"$b{-}4$", r"$b{-}3$", r"$b{-}2$", r"$b{-}1$"]
    vals = ["0", r"$b{-}1$", "2", "3", "4", "", r"$b{-}4$", r"$b{-}3$", r"$b{-}2$", "1"]
    E = 5  # ellipsis slot
    panel(ax2, labs,
          p_edges=[(0, 1), (1, 2), (8, 9)],
          v_edges=[(0, 9), (2, 9), (1, 8)],
          doubled=[(2, 3), (3, 4), (6, 7), (7, 8)],
          dashed_doubled=[(4, E), (E, 6)],
          ellipsis_slot=E, value_labels=vals)

    for ax, tag in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(0.0, 1.0, tag, transform=ax.transAxes, fontsize=9, va="top", color=INK)
    fig.tight_layout(h_pad=1.2)
    fig.savefig("fig_gadget_Tb.pdf")
    plt.close(fig)


def fig_Tpb():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 3.6))

    # (a) concrete T'_9 = (0,1,8,3,4,5,6,7,2)
    T9 = [0, 1, 8, 3, 4, 5, 6, 7, 2]
    p, v, d = concrete_edges(T9)
    panel(ax1, [str(i) for i in range(9)], p, v, d,
          value_labels=[str(x) for x in T9])

    # (b) general odd b: slots 0..4, ellipsis, b-3..b-1
    labs = ["0", "1", "2", "3", "4", "", r"$b{-}3$", r"$b{-}2$", r"$b{-}1$"]
    vals = ["0", "1", r"$b{-}1$", "3", "4", "", r"$b{-}3$", r"$b{-}2$", "2"]
    E = 5
    panel(ax2, labs,
          p_edges=[(1, 2), (2, 3), (7, 8)],
          v_edges=[(1, 8), (3, 8), (2, 7)],
          doubled=[(0, 1), (3, 4), (6, 7)],
          dashed_doubled=[(4, E), (E, 6)],
          ellipsis_slot=E, value_labels=vals)

    for ax, tag in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(0.0, 1.0, tag, transform=ax.transAxes, fontsize=9, va="top", color=INK)
    fig.tight_layout(h_pad=1.2)
    fig.savefig("fig_gadget_Tpb.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_Tb()
    fig_Tpb()
    print("wrote fig_gadget_Tb.pdf, fig_gadget_Tpb.pdf")
