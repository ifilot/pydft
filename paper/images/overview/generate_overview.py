#!/usr/bin/env python3
"""
Generate a deterministic overview figure for the PyDFT paper.

The figure is intentionally built from Matplotlib primitives rather than AI
generated imagery. The small panels are placeholders/schematics that can later
be replaced by real PyDFT output images while preserving the layout.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap


OUTDIR = Path(__file__).resolve().parent
ASSETDIR = OUTDIR / "assets"


@dataclass(frozen=True)
class Palette:
    blue: str = "#2457a6"
    blue_dark: str = "#123b7a"
    blue_light: str = "#eaf2ff"
    green: str = "#178246"
    green_dark: str = "#0e5f32"
    green_light: str = "#edf8f1"
    orange: str = "#e85d04"
    orange_dark: str = "#b94300"
    orange_light: str = "#fff3e9"
    purple: str = "#6d4c9f"
    ink: str = "#1f2933"
    muted: str = "#617085"
    line: str = "#d7dee9"
    bg: str = "#fbfcfe"


PAL = Palette()


ASSETS = {
    "molecule": ASSETDIR / "molecule.png",
    "charge": ASSETDIR / "charge.png",
    "basis": ASSETDIR / "basis.png",
    "pyqint": ASSETDIR / "pyqint.png",
    "numpy": ASSETDIR / "numpy.png",
    "pylebedev": ASSETDIR / "pylebedev.png",
}


def rounded(ax, x, y, w, h, fc="white", ec=PAL.line, lw=1.1, r=1.2, z=1):
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={r}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        zorder=z,
    )
    ax.add_patch(box)
    return box


def arrow(
    ax,
    x1,
    y1,
    x2,
    y2,
    color=PAL.blue,
    lw=1.8,
    style="-|>",
    z=5,
    connectionstyle="arc3,rad=0.0",
):
    arr = patches.FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle=style,
        mutation_scale=16,
        linewidth=lw,
        color=color,
        shrinkA=3,
        shrinkB=3,
        zorder=z,
        connectionstyle=connectionstyle,
    )
    ax.add_patch(arr)
    return arr


def elbow_arrow(ax, x1, y1, x2, y2, angle_a, angle_b, color=PAL.blue, lw=1.8, rad=8):
    return arrow(
        ax,
        x1,
        y1,
        x2,
        y2,
        color=color,
        lw=lw,
        connectionstyle=f"angle,angleA={angle_a},angleB={angle_b},rad={rad}",
    )


def badge(ax, x, y, label, color, fc="white", w=4.7, h=2.0):
    rounded(ax, x - w / 2, y - h / 2, w, h, fc=fc, ec=color, lw=1.0, r=0.45, z=8)
    text(ax, x, y, label, size=7.8, weight="bold", color=color, z=9)


def text(
    ax,
    x,
    y,
    s,
    size=9,
    weight="normal",
    color=PAL.ink,
    ha="center",
    va="center",
    z=10,
    linespacing=0.95,
):
    return ax.text(
        x,
        y,
        s,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        family="DejaVu Sans",
        linespacing=linespacing,
        zorder=z,
    )


def layer(ax, n, title, x, y, w, h, color, light):
    rounded(ax, x, y, w, h, fc=light, ec=color, lw=1.5, r=1.6, z=0)
    circ = patches.Circle((x + 3.0, y + h - 3.0), 1.7, facecolor=color, edgecolor=color, zorder=3)
    ax.add_patch(circ)
    text(ax, x + 3.0, y + h - 3.0, str(n), size=13, weight="bold", color="white")
    text(ax, x + 6.5, y + h - 3.0, title, size=15, weight="bold", color=color, ha="left")


def module(ax, x, y, w, h, title, ec=PAL.blue, fc="white", title_size=7.2):
    rounded(ax, x, y, w, h, fc=fc, ec=ec, lw=1.0, r=0.8, z=2)
    title_band = min(5.2, h * 0.38)
    text(
        ax,
        x + w / 2,
        y + h - title_band / 2,
        title,
        size=title_size,
        weight="bold",
        color=PAL.ink,
        linespacing=1.08,
    )
    return x + 0.8, y + 1.0, w - 1.6, h - title_band - 1.6


def draw_molecule(ax, x, y, w, h):
    pts = np.array(
        [
            [0.20, 0.38],
            [0.35, 0.55],
            [0.53, 0.48],
            [0.68, 0.62],
            [0.78, 0.42],
            [0.32, 0.25],
            [0.58, 0.25],
            [0.10, 0.56],
        ]
    )
    bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (1, 5), (2, 6), (0, 7)]
    for i, j in bonds:
        ax.plot(
            [x + pts[i, 0] * w, x + pts[j, 0] * w],
            [y + pts[i, 1] * h, y + pts[j, 1] * h],
            color="#525b66",
            lw=1.6,
            zorder=5,
        )
    colors = ["#f2f4f7", "#555b61", "#555b61", "#555b61", "#d63b2c", "#f2f4f7", "#f2f4f7", "#f2f4f7"]
    sizes = [45, 95, 90, 82, 95, 45, 45, 45]
    for p, c, s in zip(pts, colors, sizes):
        ax.scatter(x + p[0] * w, y + p[1] * h, s=s, c=c, edgecolors="#333", linewidths=0.8, zorder=6)


def draw_basis(ax, x, y, w, h, color=PAL.blue):
    t = np.linspace(-3.0, 3.0, 180)
    for shift, alpha in [(-1.15, 0.45), (0.0, 0.55), (1.15, 0.45)]:
        yy = np.exp(-0.8 * (t - shift) ** 2)
        xx = x + (t + 3.0) / 6.0 * w
        yy = y + 0.20 * h + yy * 0.52 * h
        ax.fill_between(xx, y + 0.20 * h, yy, color=color, alpha=alpha, zorder=4)
        ax.plot(xx, yy, color=color, lw=1.2, zorder=5)
    ax.plot([x, x + w], [y + 0.20 * h, y + 0.20 * h], color=PAL.ink, lw=1.0, zorder=5)


def draw_matrix(ax, x, y, w, h, cmap="Blues", seed=2):
    rng = np.random.default_rng(seed)
    n = 7
    mat = rng.random((n, n))
    mat = 0.5 * (mat + mat.T)
    for i in range(n):
        for j in range(n):
            ax.add_patch(
                patches.Rectangle(
                    (x + j * w / n, y + (n - 1 - i) * h / n),
                    w / n,
                    h / n,
                    facecolor=plt.get_cmap(cmap)(0.20 + 0.75 * mat[i, j]),
                    edgecolor="white",
                    linewidth=0.25,
                    zorder=4,
                )
            )


def draw_sphere_grid(ax, x, y, w, h, color=PAL.blue):
    cx, cy = x + w / 2, y + h / 2
    r = min(w, h) * 0.42
    ax.add_patch(patches.Circle((cx, cy), r, facecolor="#f5f9ff", edgecolor=color, lw=1.1, zorder=3))
    for k in [-0.65, -0.32, 0.0, 0.32, 0.65]:
        ax.add_patch(patches.Ellipse((cx, cy), 2 * r, 2 * r * np.sqrt(1 - k * k), fill=False, edgecolor=color, lw=0.45, alpha=0.65, zorder=4))
    for angle in np.linspace(0, np.pi, 6, endpoint=False):
        ax.add_patch(patches.Ellipse((cx, cy), 2 * r * np.cos(angle) ** 2 + 0.1, 2 * r, angle=np.degrees(angle), fill=False, edgecolor=color, lw=0.4, alpha=0.45, zorder=4))
    for theta in np.linspace(0, 2 * np.pi, 14, endpoint=False):
        for rr in [0.35, 0.68, 0.92]:
            ax.scatter(cx + r * rr * np.cos(theta), cy + r * rr * np.sin(theta), s=7, c=color, zorder=5)


def draw_grid_dots(ax, x, y, w, h, color=PAL.green):
    draw_molecule(ax, x + 0.18 * w, y + 0.20 * h, 0.58 * w, 0.58 * h)
    cx, cy = x + 0.5 * w, y + 0.52 * h
    for rr in np.linspace(0.15, 0.46, 5):
        for th in np.linspace(0, 2 * np.pi, 16, endpoint=False):
            ax.scatter(cx + rr * w * np.cos(th), cy + rr * h * np.sin(th), s=3.8, c=color, zorder=3)


def draw_fuzzy(ax, x, y, w, h):
    draw_molecule(ax, x + 0.18 * w, y + 0.20 * h, 0.58 * w, 0.58 * h)
    centers = [(0.28, 0.45), (0.50, 0.55), (0.70, 0.46)]
    colors = [PAL.green, "#8ccf9f", "#d7a246"]
    for (cx, cy), c in zip(centers, colors):
        for rr, alpha in [(0.35, 0.10), (0.26, 0.14), (0.16, 0.18)]:
            ax.add_patch(
                patches.Ellipse(
                    (x + cx * w, y + cy * h),
                    rr * w,
                    rr * h,
                    facecolor=c,
                    edgecolor="none",
                    alpha=alpha,
                    zorder=2,
                )
            )
    for th in np.linspace(0, 2 * np.pi, 7, endpoint=False):
        ax.plot(
            [x + 0.5 * w, x + 0.5 * w + 0.42 * w * np.cos(th)],
            [y + 0.5 * h, y + 0.5 * h + 0.42 * h * np.sin(th)],
            color="#78828f",
            lw=0.45,
            ls="--",
            zorder=3,
        )


def draw_curve(ax, x, y, w, h, color=PAL.green):
    ax.plot([x + 0.10 * w, x + 0.10 * w], [y + 0.17 * h, y + 0.86 * h], color=PAL.ink, lw=0.8)
    ax.plot([x + 0.10 * w, x + 0.92 * w], [y + 0.17 * h, y + 0.17 * h], color=PAL.ink, lw=0.8)
    t = np.linspace(0, 1, 9)
    xx = x + (0.14 + 0.74 * t) * w
    yy = y + (0.18 + 0.55 * np.sin(np.pi * t) ** 0.7) * h
    ax.plot(xx, yy, color=color, lw=1.2)
    ax.scatter(xx, yy, s=18, c="#f0f7ff", edgecolors=color, linewidths=1.0, zorder=6)


def draw_density(ax, x, y, w, h, color="#279989"):
    nx = ny = 90
    xx, yy = np.meshgrid(np.linspace(-2.0, 2.0, nx), np.linspace(-2.0, 2.0, ny))
    z = np.exp(-((xx - 0.55) ** 2 + (yy + 0.1) ** 2) * 1.8) + 0.8 * np.exp(-((xx + 0.65) ** 2 + (yy - 0.2) ** 2) * 1.4)
    levels = np.linspace(z.min() + 0.08, z.max(), 8)
    ax.contourf(x + (xx + 2) / 4 * w, y + (yy + 2) / 4 * h, z, levels=levels, cmap="GnBu", alpha=0.75, zorder=3)
    ax.contour(x + (xx + 2) / 4 * w, y + (yy + 2) / 4 * h, z, levels=levels, colors="#27606a", linewidths=0.45, alpha=0.75, zorder=4)
    ax.arrow(x + 0.52 * w, y + 0.52 * h, 0.16 * w, 0.10 * h, head_width=0.22, head_length=0.22, color=color, lw=0.8, zorder=5)


def draw_hartree(ax, x, y, w, h):
    nx = ny = 80
    xx, yy = np.meshgrid(np.linspace(-2.0, 2.0, nx), np.linspace(-2.0, 2.0, ny))
    z = np.exp(-(xx**2 + (yy - 0.55) ** 2) * 1.3) + 0.7 * np.exp(-((xx * 1.2) ** 2 + (yy + 0.9) ** 2) * 1.8)
    ax.contourf(x + (xx + 2) / 4 * w, y + (yy + 2) / 4 * h, z, levels=10, cmap="Purples", alpha=0.85, zorder=3)
    ax.contour(x + (xx + 2) / 4 * w, y + (yy + 2) / 4 * h, z, levels=8, colors="#4a3a75", linewidths=0.35, zorder=4)


def draw_xc(ax, x, y, w, h):
    ax.plot([x + 0.12 * w, x + 0.12 * w], [y + 0.20 * h, y + 0.84 * h], color=PAL.ink, lw=0.8)
    ax.plot([x + 0.12 * w, x + 0.90 * w], [y + 0.20 * h, y + 0.20 * h], color=PAL.ink, lw=0.8)
    t = np.linspace(0.08, 1, 160)
    y1 = 0.72 - 0.38 * np.sqrt(t)
    y2 = 0.60 - 0.28 * t**0.75 - 0.07 * np.sin(4 * t)
    ax.plot(x + (0.12 + 0.75 * t) * w, y + y1 * h, color=PAL.orange, lw=1.4)
    ax.plot(x + (0.12 + 0.75 * t) * w, y + y2 * h, color=PAL.blue, lw=1.2, ls="--")
    text(ax, x + 0.70 * w, y + 0.72 * h, "LDA", size=6.5, color=PAL.orange)
    text(ax, x + 0.70 * w, y + 0.55 * h, "PBE", size=6.5, color=PAL.blue)


def draw_orbital(ax, x, y, w, h):
    lobes = [
        (0.34, 0.58, 0.20, 0.34, "#ce3f32"),
        (0.66, 0.42, 0.20, 0.34, "#ce3f32"),
        (0.34, 0.42, 0.18, 0.27, PAL.blue),
        (0.66, 0.58, 0.18, 0.27, PAL.blue),
    ]
    for cx, cy, ww, hh, c in lobes:
        ax.add_patch(patches.Ellipse((x + cx * w, y + cy * h), ww * w, hh * h, angle=25, facecolor=c, edgecolor="#333", lw=0.4, alpha=0.9, zorder=4))
    ax.scatter([x + 0.48 * w, x + 0.54 * w], [y + 0.50 * h, y + 0.50 * h], s=35, c="#555", zorder=5)


def draw_pie(ax, x, y, w, h):
    center = (x + 0.5 * w, y + 0.50 * h)
    vals = np.array([0.26, 0.24, 0.18, 0.17, 0.15])
    colors = [PAL.blue, PAL.orange, PAL.green, "#8b6bb8", "#8895a7"]
    theta = 90
    r = min(w, h) * 0.34
    for v, c in zip(vals, colors):
        ax.add_patch(patches.Wedge(center, r, theta, theta + 360 * v, facecolor=c, edgecolor="white", lw=0.8, zorder=4))
        theta += 360 * v


def draw_sliders(ax, x, y, w, h):
    labels = ["basis", "grid", "XC", "tol", "charge"]
    for i, lab in enumerate(labels):
        yy = y + h * (0.82 - i * 0.16)
        text(ax, x + 0.18 * w, yy, lab, size=6.8, color=PAL.ink, ha="left")
        ax.plot([x + 0.48 * w, x + 0.90 * w], [yy, yy], color="#8995a5", lw=1.0, zorder=4)
        ax.scatter(x + (0.56 + 0.08 * (i % 3)) * w, yy, s=28, c=PAL.blue, edgecolors=PAL.blue_dark, linewidths=0.5, zorder=5)


def draw_mini(ax, kind, x, y, w, h):
    if kind == "molecule":
        draw_molecule(ax, x, y, w, h)
    elif kind == "basis":
        draw_basis(ax, x, y, w, h)
    elif kind == "matrix":
        draw_matrix(ax, x + 0.12 * w, y + 0.12 * h, 0.76 * w, 0.64 * h)
    elif kind == "sphere":
        draw_sphere_grid(ax, x, y, w, h)
    elif kind == "grid":
        draw_grid_dots(ax, x, y, w, h)
    elif kind == "fuzzy":
        draw_fuzzy(ax, x, y, w, h)
    elif kind == "curve":
        draw_curve(ax, x, y, w, h)
    elif kind == "density":
        draw_density(ax, x + 0.08 * w, y + 0.08 * h, 0.84 * w, 0.68 * h)
    elif kind == "hartree":
        draw_hartree(ax, x + 0.10 * w, y + 0.10 * h, 0.80 * w, 0.70 * h)
    elif kind == "xc":
        draw_xc(ax, x, y, w, h)
    elif kind == "orbital":
        draw_orbital(ax, x, y, w, h)
    elif kind == "pie":
        draw_pie(ax, x, y, w, h)
    elif kind == "sliders":
        draw_sliders(ax, x, y, w, h)


def draw_asset(ax, name, x, y, w, h):
    img = mpimg.imread(ASSETS[name])
    ax.imshow(img, extent=(x, x + w, y, y + h), aspect="auto", zorder=4)


def build():
    fig, ax = plt.subplots(figsize=(12, 16.5), dpi=180)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 140)
    ax.axis("off")
    fig.patch.set_facecolor(PAL.bg)
    ax.set_facecolor(PAL.bg)

    layer(ax, 1, "Inputs and dependencies", 1.2, 112.0, 97.6, 26.0, PAL.blue, PAL.blue_light)
    layer(ax, 2, "Core PyDFT construction", 1.2, 82.0, 97.6, 26.0, PAL.green, PAL.green_light)
    layer(ax, 3, "Kohn-Sham iteration", 1.2, 42.5, 97.6, 35.2, PAL.blue, PAL.blue_light)
    layer(ax, 4, "Exposed learning outputs", 1.2, 8.5, 97.6, 29.0, PAL.orange, PAL.orange_light)

    # Layer 1
    left_x, left_y, left_w, left_h = 3.2, 115.0, 47.0, 17.4
    right_x, right_y, right_w, right_h = 54.0, 115.0, 42.8, 17.4
    rounded(ax, left_x, left_y, left_w, left_h, fc="white", ec=PAL.blue, lw=1.0, r=1.0)
    rounded(ax, right_x, right_y, right_w, right_h, fc="white", ec=PAL.blue, lw=1.0, r=1.0)
    for i, (title, kind, sub) in enumerate(
        [
            ("Molecule", "molecule", "geometry"),
            ("Charge", "charge", "closed shell"),
            ("Basis set", "basis", "STO-3G, 6-31G(d)"),
        ]
    ):
        x = left_x + 2 + i * 15.0
        text(ax, x + 6.6, left_y + left_h - 2.3, title, size=8.5, weight="bold")
        draw_asset(ax, kind, x + 1.0, left_y + 4.35, 11.2, 8.4)
        text(ax, x + 6.6, left_y + 2.1, sub, size=7.0, color=PAL.ink)

    for i, (title, kind, sub) in enumerate(
        [
            ("PyQInt", "pyqint", "analytic integrals"),
            ("NumPy", "numpy", "arrays"),
            ("PyLebedev", "pylebedev", "angular rules"),
        ]
    ):
        x = right_x + 2.0 + i * 13.4
        text(ax, x + 5.8, right_y + right_h - 2.3, title, size=8.5, weight="bold", color=PAL.blue_dark)
        draw_asset(ax, kind, x + 0.7, right_y + 4.35, 10.2, 8.4)
        text(ax, x + 5.8, right_y + 2.2, sub, size=6.8, color=PAL.ink)

    arrow(ax, 50, 111.7, 50, 108.7, color=PAL.blue, lw=2.0)

    # Layer 2
    core_items = [
        ("Basis\namplitudes", "basis"),
        ("$S_{\\mu\\nu}, T_{\\mu\\nu}, V_{\\mu\\nu}$", "matrix"),
        ("Becke grid\n$\\{\\mathbf{r}_i,w_i\\}$", "grid"),
        ("Fuzzy weights\n$w_A(\\mathbf{r})$", "fuzzy"),
        ("Lebedev\n$\\int d\\Omega$", "sphere"),
        ("Radial grid\n$r_k$", "curve"),
    ]
    x0, y0, gap = 3.2, 85.1, 1.7
    mw = (93.6 - gap * 5) / 6
    for i, (title, kind) in enumerate(core_items):
        x = x0 + i * (mw + gap)
        cx, cy, cw, ch = module(ax, x, y0, mw, 16.8, title, ec=PAL.green, title_size=6.8)
        draw_mini(ax, kind, cx, cy + 0.2, cw, ch - 0.2)

    arrow(ax, 50, 81.9, 50, 78.0, color=PAL.green, lw=2.0)

    # Layer 3
    iter_items = [
        ("Density matrix\n$P_{\\mu\\nu}$", "matrix"),
        ("$\\rho(\\mathbf{r})$\n$\\nabla\\rho(\\mathbf{r})$", "density"),
        ("Hartree\n$v_H(\\mathbf{r})$", "hartree"),
        ("XC\n$v_{xc}^{\\mathrm{LDA/PBE}}$", "xc"),
        ("Build Fock\n$F=H+J+V_{xc}$", "matrix"),
        ("Solve KS\n$FC=SC\\epsilon$", None),
        ("Update density\n$P=2CC^T$", "matrix"),
    ]
    ix0, iy0, igap = 3.2, 54.0, 2.2
    iw = (93.6 - igap * 6) / 7
    centers = []
    cards = []
    for i, (title, kind) in enumerate(iter_items):
        x = ix0 + i * (iw + igap)
        cards.append((x, iy0, iw, 16.2))
        cx, cy, cw, ch = module(ax, x, iy0, iw, 16.2, title, ec=PAL.blue, title_size=6.5)
        if kind:
            draw_mini(ax, kind, cx, cy + 0.2, cw, ch - 0.2)
        else:
            for k in range(3):
                yy = cy + 1.8 + k * 2.2
                ax.plot([x + 2.4, x + iw - 2.4], [yy, yy], color=PAL.ink, lw=0.8)
                ax.scatter([x + iw / 2], [yy], s=20 + 10 * k, c=PAL.blue, zorder=5)
        centers.append((x + iw / 2, iy0 + 8.0))
        if i > 0:
            prev_x, prev_y, prev_w, prev_h = cards[i - 1]
            curr_x, curr_y, curr_w, curr_h = cards[i]
            gap_mid_y = curr_y + curr_h / 2
            arrow(
                ax,
                prev_x + prev_w + 0.10,
                gap_mid_y,
                curr_x - 0.10,
                gap_mid_y,
                color=PAL.blue,
                lw=1.5,
            )

    solve_x, solve_y, solve_w, _solve_h = cards[5]
    diamond_x = solve_x + solve_w / 2
    diamond_y = 46.2
    diamond_w = 8.2
    diamond_h = 6.8
    rounded(
        ax,
        diamond_x - diamond_w / 2,
        diamond_y - diamond_h / 2,
        diamond_w,
        diamond_h,
        fc="white",
        ec=PAL.blue,
        lw=1.0,
        r=0.25,
        z=4,
    )
    text(ax, diamond_x, diamond_y + 0.9, "Converged?", size=7.2, weight="bold")
    text(ax, diamond_x, diamond_y - 1.0, "$|\\Delta E|, |\\Delta P| < \\tau$", size=6.4, color=PAL.muted)
    first_x, first_y, first_w, _first_h = cards[0]
    arrow(
        ax,
        diamond_x,
        solve_y - 0.15,
        diamond_x,
        diamond_y + diamond_h / 2,
        color=PAL.blue,
        lw=1.5,
    )
    elbow_arrow(
        ax,
        diamond_x - diamond_w / 2,
        diamond_y - 0.1,
        first_x + first_w / 2,
        first_y - 0.1,
        angle_a=180,
        angle_b=-90,
        color=PAL.blue,
        lw=1.5,
        rad=9,
    )
    arrow(ax, diamond_x, diamond_y - diamond_h / 2, diamond_x, 35.8, color=PAL.green, lw=2.0)
    badge(ax, 35.0, diamond_y - 0.1, "No", "#c52020", fc="#fff5f5", w=4.8, h=2.0)
    badge(ax, diamond_x, 39.4, "Yes", PAL.green, fc="#f0faf4", w=5.2, h=2.0)

    # Layer 4
    output_items = [
        ("Matrices\n$F,P,S$", "matrix"),
        ("Energy terms\n$E_i$", "pie"),
        ("Orbitals\n$\\psi_i$", "orbital"),
        ("Density maps\n$\\rho(\\mathbf{r})$", "density"),
        ("Becke cells\n$w_A(\\mathbf{r})$", "fuzzy"),
        ("Spherical coeff.\n$\\rho_{klm}$", "curve"),
        ("Timing\nbreakdown", "matrix"),
        ("Interactive\nexploration", "sliders"),
    ]
    ox0, oy0, ogap = 3.2, 13.2, 1.15
    ow = (93.6 - ogap * 7) / 8
    oh = 18.2
    for i, (title, kind) in enumerate(output_items):
        x = ox0 + i * (ow + ogap)
        cx, cy, cw, ch = module(
            ax,
            x,
            oy0,
            ow,
            oh,
            title,
            ec=PAL.orange,
            fc="white",
            title_size=6.35,
        )
        draw_mini(ax, kind, cx, cy + 0.3, cw, ch - 0.3)

    # Figure legend.
    yleg = 3.6
    arrow(ax, 5.0, yleg, 9.0, yleg, color=PAL.blue, lw=1.8)
    text(ax, 16.0, yleg, "data / control flow", size=7.2, color=PAL.ink)
    arrow(ax, 31.0, yleg, 35.0, yleg, color=PAL.green, lw=1.8)
    text(ax, 42.5, yleg, "accepted workflow", size=7.2, color=PAL.ink)
    rounded(ax, 55.5, yleg - 1.0, 3.7, 2.0, fc="white", ec=PAL.green, lw=1.0, r=0.3)
    text(ax, 66.0, yleg, "construction step", size=7.2, color=PAL.ink)
    rounded(ax, 78.0, yleg - 1.0, 3.7, 2.0, fc="white", ec=PAL.orange, lw=1.0, r=0.3)
    text(ax, 88.2, yleg, "learning output", size=7.2, color=PAL.ink)

    png = OUTDIR / "overview-draft.png"
    pdf = OUTDIR / "overview-draft.pdf"
    svg = OUTDIR / "overview-draft.svg"
    fig.savefig(png, dpi=240, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(svg, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"Wrote {png}")
    print(f"Wrote {pdf}")
    print(f"Wrote {svg}")


if __name__ == "__main__":
    build()
