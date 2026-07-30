#!/usr/bin/env python3
"""Generate the gadget figures (E12) in the style of the union_viz.html teaching tool
(misc/PC_Ham_Teaching/tools/union_viz.html in the Periapsis repo).

Conventions (matching the tool and the paper's terminology):
  - vertices = positions 0..b-1 as dark dots on a gray baseline;
  - BLUE arcs above the line = P-edges {i,i+1};
  - RED arcs below the line = V-edges, arc width proportional to the position jump;
  - PURPLE paired arcs (one above, one below) = doubled pairs;
  - labels in the tool's monospace font (ui-monospace/Consolas), position labels bold;
  - general-b panels elide the middle of the doubled run with a centered ellipsis;
    SOLID purple stubs (~1/3 of the gap each side) indicate the run continues,
    leaving the ellipsis visible.

Outputs: fig_gadget_Tb.pdf  (T_8 concrete + general even b)
         fig_gadget_Tpb.pdf (T'_9 concrete + general odd b)
Run from this directory:  python gen_gadget_figs.py
"""

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "monospace",
    "font.monospace": ["Consolas", "Cascadia Mono", "DejaVu Sans Mono"],
    "pdf.fonttype": 42,  # embed as TrueType, not Type 3
})
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

BLUE, RED, DOUBLE, INK, MUTED, LINE = "#2196f3", "#e53935", "#8e24aa", "#1a1a1a", "#666666", "#333333"
VSCALE = 0.30   # arch peak = VSCALE * span (tool default 0.42; reduced for print)
MINUS = "−"
ELLIPSIS = "···"  # middle dots; present in Consolas
STUB_T = 0.40   # Bezier parameter cut: horizontal extent 3t^2-2t^3 ~ 0.35 of the gap


def arch_pts(xa, xb, y0, direction, vscale):
    h = direction * (4.0 / 3.0) * vscale * abs(xb - xa)
    return [(xa, y0), (xa, y0 + h), (xb, y0 + h), (xb, y0)]


def _lerp(p, q, t):
    return (p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t)


def _split(pts, t):
    """De Casteljau: control points of the [0,t] and [t,1] halves of a cubic."""
    p0, p1, p2, p3 = pts
    a, b, c = _lerp(p0, p1, t), _lerp(p1, p2, t), _lerp(p2, p3, t)
    d, e = _lerp(a, b, t), _lerp(b, c, t)
    f = _lerp(d, e, t)
    return (p0, a, d, f), (f, e, c, p3)


def _draw_cubic(ax, pts, color, lw):
    path = Path(list(pts), [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, fill=False, edgecolor=color, lw=lw,
                           capstyle="round", alpha=0.95))


def arch(ax, xa, xb, y0, direction, color, lw, vscale=VSCALE):
    """Cubic-Bezier arch from (xa,y0) to (xb,y0); direction +1 up, -1 down."""
    _draw_cubic(ax, arch_pts(xa, xb, y0, direction, vscale), color, lw)


def stub(ax, xa, xb, y0, direction, color, lw, at_start, vscale=VSCALE):
    """Truncated arch: the first (at_start) or last ~1/3 of the arch xa->xb."""
    pts = arch_pts(xa, xb, y0, direction, vscale)
    left, right = _split(pts, STUB_T if at_start else 1.0 - STUB_T)
    _draw_cubic(ax, left if at_start else right, color, lw)


def panel(ax, slots, p_edges, v_edges, doubled, stub_doubled=(), ellipsis_slot=None,
          value_labels=None, vscale=VSCALE):
    """Draw one H_sigma panel.

    slots: list of position-label strings ('' allowed), one per x-slot.
    p_edges / v_edges: (i, j) slot pairs for P-only (top, blue) / V-only (bottom, red).
    doubled: (i, j) slot pairs drawn as a purple arc pair (top + bottom).
    stub_doubled: like doubled, but one endpoint is the ellipsis slot: only the ~1/3
      of the arc pair nearest the real vertex is drawn (solid), so the run visibly
      continues into the elision without covering the ellipsis.
    ellipsis_slot: slot index that holds a centered ellipsis instead of a vertex.
    value_labels: optional list of value-label strings under each vertex.
    """
    n = len(slots)
    top = max((abs(b - a) for a, b in list(p_edges) + list(doubled) + list(stub_doubled)), default=1)
    bot = max((abs(b - a) for a, b in list(v_edges) + list(doubled) + list(stub_doubled)), default=1)
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
    for a, b in stub_doubled:
        at_start = (b == ellipsis_slot)  # stub leaves the real vertex toward the elision
        for direction in (+1, -1):
            stub(ax, a, b, 0, direction, DOUBLE, 2.0, at_start)

    for i, lab in enumerate(slots):
        if i == ellipsis_slot:
            ax.text(i, 0, ELLIPSIS, ha="center", va="center", fontsize=10,
                    color=INK, zorder=3)
            continue
        ax.plot([i], [0], "o", ms=3.6, color=LINE, zorder=3)
        ax.text(i, vscale * top + 0.34, lab, ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=INK, zorder=3)
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


def bm(k):
    """Label string 'b-k' with a real minus sign."""
    return "b" + MINUS + str(k)


def fig_Tb():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.0, 3.6))

    # (a) concrete T_8 = (0,7,2,3,4,5,6,1)
    T8 = [0, 7, 2, 3, 4, 5, 6, 1]
    p, v, d = concrete_edges(T8)
    panel(ax1, [str(i) for i in range(8)], p, v, d,
          value_labels=[str(x) for x in T8])

    # (b) general even b: slots 0..4, ellipsis, b-4..b-1
    labs = ["0", "1", "2", "3", "4", "", bm(4), bm(3), bm(2), bm(1)]
    vals = ["0", bm(1), "2", "3", "4", "", bm(4), bm(3), bm(2), "1"]
    E = 5  # ellipsis slot
    panel(ax2, labs,
          p_edges=[(0, 1), (1, 2), (8, 9)],
          v_edges=[(0, 9), (2, 9), (1, 8)],
          doubled=[(2, 3), (3, 4), (6, 7), (7, 8)],
          stub_doubled=[(4, E), (E, 6)],
          ellipsis_slot=E, value_labels=vals)

    for ax, tag in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(-0.03, 1.0, tag, transform=ax.transAxes, fontsize=9, va="top",
                ha="left", color=INK)
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
    labs = ["0", "1", "2", "3", "4", "", bm(3), bm(2), bm(1)]
    vals = ["0", "1", bm(1), "3", "4", "", bm(3), bm(2), "2"]
    E = 5
    panel(ax2, labs,
          p_edges=[(1, 2), (2, 3), (7, 8)],
          v_edges=[(1, 8), (3, 8), (2, 7)],
          doubled=[(0, 1), (3, 4), (6, 7)],
          stub_doubled=[(4, E), (E, 6)],
          ellipsis_slot=E, value_labels=vals)

    for ax, tag in ((ax1, "(a)"), (ax2, "(b)")):
        ax.text(-0.03, 1.0, tag, transform=ax.transAxes, fontsize=9, va="top",
                ha="left", color=INK)
    fig.tight_layout(h_pad=1.2)
    fig.savefig("fig_gadget_Tpb.pdf")
    plt.close(fig)


if __name__ == "__main__":
    fig_Tb()
    fig_Tpb()
    print("wrote fig_gadget_Tb.pdf, fig_gadget_Tpb.pdf")
