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
import re

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from matplotlib import patches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MplPath
from matplotlib.transforms import Affine2D


OUTDIR = Path(__file__).resolve().parent
ASSETDIR = OUTDIR / "assets"


@dataclass(frozen=True)
class Palette:
    blue: str = "#2457a6"
    blue_dark: str = "#123b7a"
    blue_light: str = "#eaf2ff"
    azure: str = "#1282c4"
    azure_dark: str = "#075985"
    azure_light: str = "#e8f7ff"
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
    "basis-amplitudes": ASSETDIR / "basis-amplitudes.png",
    "one-electron-integrals": ASSETDIR / "one-electron-integrals.png",
    "becke-grid": ASSETDIR / "becke-grid.png",
    "fuzzy-weights": ASSETDIR / "fuzzy-weights.png",
    "lebedev-grid": ASSETDIR / "lebedev-grid.png",
    "radial-grid": ASSETDIR / "radial-grid.png",
    "density-matrix": ASSETDIR / "density-matrix.png",
    "density-gradient": ASSETDIR / "density-gradient.png",
    "hartree-potential": ASSETDIR / "hartree-potential.png",
    "xc-potential": ASSETDIR / "xc-potential.png",
    "fock-matrix-build": ASSETDIR / "fock-matrix-build.png",
    "solve-ks": ASSETDIR / "solve-ks.png",
    "update-density": ASSETDIR / "update-density.png",
    "output-matrices": ASSETDIR / "output-matrices.png",
    "energy-terms": ASSETDIR / "energy-terms.png",
    "orbitals": ASSETDIR / "orbitals.png",
    "density-maps": ASSETDIR / "density-maps.png",
    "becke-cells": ASSETDIR / "becke-cells.png",
    "spherical-coefficients": ASSETDIR / "spherical-coefficients.png",
    "timing-breakdown": ASSETDIR / "timing-breakdown.png",
    "interactive-exploration": ASSETDIR / "interactive-exploration.png",
}


# Font Awesome Free 7.2.0 SVG path data.
FA_ICONS = {
    "inbox": {
        "viewbox": (512, 512),
        "path": "M91.8 32C59.9 32 32.9 55.4 28.4 86.9L.6 281.2c-.4 3-.6 6-.6 9.1L0 416c0 35.3 28.7 64 64 64l384 0c35.3 0 64-28.7 64-64l0-125.7c0-3-.2-6.1-.6-9.1L483.6 86.9C479.1 55.4 452.1 32 420.2 32L91.8 32zm0 64l328.5 0 27.4 192-59.9 0c-12.1 0-23.2 6.8-28.6 17.7l-14.3 28.6c-5.4 10.8-16.5 17.7-28.6 17.7l-120.4 0c-12.1 0-23.2-6.8-28.6-17.7l-14.3-28.6c-5.4-10.8-16.5-17.7-28.6-17.7L64.3 288 91.8 96z",
    },
    "gear": {
        "viewbox": (512, 512),
        "path": "M195.1 9.5C198.1-5.3 211.2-16 226.4-16l59.8 0c15.2 0 28.3 10.7 31.3 25.5L332 79.5c14.1 6 27.3 13.7 39.3 22.8l67.8-22.5c14.4-4.8 30.2 1.2 37.8 14.4l29.9 51.8c7.6 13.2 4.9 29.8-6.5 39.9L447 233.3c.9 7.4 1.3 15 1.3 22.7s-.5 15.3-1.3 22.7l53.4 47.5c11.4 10.1 14 26.8 6.5 39.9l-29.9 51.8c-7.6 13.1-23.4 19.2-37.8 14.4l-67.8-22.5c-12.1 9.1-25.3 16.7-39.3 22.8l-14.4 69.9c-3.1 14.9-16.2 25.5-31.3 25.5l-59.8 0c-15.2 0-28.3-10.7-31.3-25.5l-14.4-69.9c-14.1-6-27.2-13.7-39.3-22.8L73.5 432.3c-14.4 4.8-30.2-1.2-37.8-14.4L5.8 366.1c-7.6-13.2-4.9-29.8 6.5-39.9l53.4-47.5c-.9-7.4-1.3-15-1.3-22.7s.5-15.3 1.3-22.7L12.3 185.8c-11.4-10.1-14-26.8-6.5-39.9L35.7 94.1c7.6-13.2 23.4-19.2 37.8-14.4l67.8 22.5c12.1-9.1 25.3-16.7 39.3-22.8L195.1 9.5zM256.3 336a80 80 0 1 0-.6-160 80 80 0 1 0 .6 160z",
    },
    "arrows-rotate": {
        "viewbox": (512, 512),
        "path": "M65.9 228.5c13.3-93 93.4-164.5 190.1-164.5 53 0 101 21.5 135.8 56.2.2.2.4.4.6.6l7.6 7.2-47.9 0c-17.7 0-32 14.3-32 32s14.3 32 32 32l128 0c17.7 0 32-14.3 32-32l0-128c0-17.7-14.3-32-32-32s-32 14.3-32 32l0 53.4-11.3-10.7C390.5 28.6 326.5 0 256 0 127 0 20.3 95.4 2.6 219.5.1 237 12.2 253.2 29.7 255.7s33.7-9.7 36.2-27.1zm443.5 64c2.5-17.5-9.7-33.7-27.1-36.2s-33.7 9.7-36.2 27.1c-13.3 93-93.4 164.5-190.1 164.5-53 0-101-21.5-135.8-56.2-.2-.2-.4-.4-.6-.6l-7.6-7.2 47.9 0c17.7 0 32-14.3 32-32s-14.3-32-32-32L32 320c-8.5 0-16.7 3.4-22.7 9.5S-.1 343.7 0 352.3l1 127c.1 17.7 14.6 31.9 32.3 31.7S65.2 496.4 65 478.7l-.4-51.5 10.7 10.1c46.3 46.1 110.2 74.7 180.7 74.7 129 0 235.7-95.4 253.4-219.5z",
    },
    "chart-line": {
        "viewbox": (512, 512),
        "path": "M64 64c0-17.7-14.3-32-32-32S0 46.3 0 64L0 400c0 44.2 35.8 80 80 80l400 0c17.7 0 32-14.3 32-32s-14.3-32-32-32L80 416c-8.8 0-16-7.2-16-16L64 64zm406.6 86.6c12.5-12.5 12.5-32.8 0-45.3s-32.8-12.5-45.3 0L320 210.7 262.6 153.4c-12.5-12.5-32.8-12.5-45.3 0l-96 96c-12.5 12.5-12.5 32.8 0 45.3s32.8 12.5 45.3 0l73.4-73.4 57.4 57.4c12.5 12.5 32.8 12.5 45.3 0l128-128z",
    },
}


def rounded(ax, x, y, w, h, fc="white", ec=PAL.line, lw=1.1, r=1.2, z=1, alpha=1.0):
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={r}",
        facecolor=fc,
        edgecolor=ec,
        linewidth=lw,
        alpha=alpha,
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


def _svg_tokens(path_data):
    token_re = r"[AaCcHhLlMmQqSsVvZz]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?"
    return re.findall(token_re, path_data)


def _is_cmd(token):
    return len(token) == 1 and token.isalpha()


def _arc_points(p0, rx, ry, angle, large_arc, sweep, p1):
    if rx == 0 or ry == 0:
        return [p1]
    x0, y0 = p0
    x1, y1 = p1
    rx, ry = abs(rx), abs(ry)
    phi = np.deg2rad(angle)
    cos_phi, sin_phi = np.cos(phi), np.sin(phi)
    dx2, dy2 = (x0 - x1) / 2.0, (y0 - y1) / 2.0
    x1p = cos_phi * dx2 + sin_phi * dy2
    y1p = -sin_phi * dx2 + cos_phi * dy2

    lam = (x1p**2) / (rx**2) + (y1p**2) / (ry**2)
    if lam > 1:
        scale = np.sqrt(lam)
        rx *= scale
        ry *= scale

    denom = rx**2 * y1p**2 + ry**2 * x1p**2
    if denom == 0:
        return [p1]
    numer = max(0.0, rx**2 * ry**2 - denom)
    coef = (-1 if large_arc == sweep else 1) * np.sqrt(numer / denom)
    cxp = coef * (rx * y1p / ry)
    cyp = coef * (-ry * x1p / rx)
    cx = cos_phi * cxp - sin_phi * cyp + (x0 + x1) / 2.0
    cy = sin_phi * cxp + cos_phi * cyp + (y0 + y1) / 2.0

    def angle_between(u, v):
        dot = np.dot(u, v)
        det = u[0] * v[1] - u[1] * v[0]
        return np.arctan2(det, dot)

    u = np.array([(x1p - cxp) / rx, (y1p - cyp) / ry])
    v = np.array([(-x1p - cxp) / rx, (-y1p - cyp) / ry])
    theta1 = angle_between(np.array([1.0, 0.0]), u)
    delta = angle_between(u, v)
    if not sweep and delta > 0:
        delta -= 2 * np.pi
    elif sweep and delta < 0:
        delta += 2 * np.pi

    n = max(4, int(np.ceil(abs(delta) / (np.pi / 12))))
    pts = []
    for i in range(1, n + 1):
        theta = theta1 + delta * i / n
        x = cx + rx * cos_phi * np.cos(theta) - ry * sin_phi * np.sin(theta)
        y = cy + rx * sin_phi * np.cos(theta) + ry * cos_phi * np.sin(theta)
        pts.append((x, y))
    return pts


def svg_path_to_mpl(path_data):
    tokens = _svg_tokens(path_data)
    verts = []
    codes = []
    i = 0
    cmd = None
    cur = np.array([0.0, 0.0])
    start = np.array([0.0, 0.0])
    last_cubic_ctrl = None

    def read_float():
        nonlocal i
        val = float(tokens[i])
        i += 1
        return val

    def has_number():
        return i < len(tokens) and not _is_cmd(tokens[i])

    while i < len(tokens):
        if _is_cmd(tokens[i]):
            cmd = tokens[i]
            i += 1
        if cmd is None:
            break

        lower = cmd.lower()
        rel = cmd.islower()
        if lower == "m":
            first = True
            while has_number():
                p = np.array([read_float(), read_float()])
                if rel:
                    p = cur + p
                codes.append(MplPath.MOVETO if first else MplPath.LINETO)
                verts.append(tuple(p))
                cur = p
                if first:
                    start = p.copy()
                    first = False
            last_cubic_ctrl = None
        elif lower == "l":
            while has_number():
                p = np.array([read_float(), read_float()])
                if rel:
                    p = cur + p
                codes.append(MplPath.LINETO)
                verts.append(tuple(p))
                cur = p
            last_cubic_ctrl = None
        elif lower == "h":
            while has_number():
                x = read_float()
                p = np.array([cur[0] + x if rel else x, cur[1]])
                codes.append(MplPath.LINETO)
                verts.append(tuple(p))
                cur = p
            last_cubic_ctrl = None
        elif lower == "v":
            while has_number():
                y = read_float()
                p = np.array([cur[0], cur[1] + y if rel else y])
                codes.append(MplPath.LINETO)
                verts.append(tuple(p))
                cur = p
            last_cubic_ctrl = None
        elif lower == "c":
            while has_number():
                pts = [np.array([read_float(), read_float()]) for _ in range(3)]
                if rel:
                    pts = [cur + p for p in pts]
                verts.extend([tuple(p) for p in pts])
                codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
                cur = pts[2]
                last_cubic_ctrl = pts[1]
        elif lower == "s":
            while has_number():
                c1 = cur if last_cubic_ctrl is None else cur + (cur - last_cubic_ctrl)
                c2 = np.array([read_float(), read_float()])
                p = np.array([read_float(), read_float()])
                if rel:
                    c2 = cur + c2
                    p = cur + p
                verts.extend([tuple(c1), tuple(c2), tuple(p)])
                codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
                cur = p
                last_cubic_ctrl = c2
        elif lower == "a":
            while has_number():
                rx, ry = read_float(), read_float()
                angle = read_float()
                large_arc, sweep = int(read_float()), int(read_float())
                p = np.array([read_float(), read_float()])
                if rel:
                    p = cur + p
                for pt in _arc_points(tuple(cur), rx, ry, angle, large_arc, sweep, tuple(p)):
                    codes.append(MplPath.LINETO)
                    verts.append(pt)
                cur = p
            last_cubic_ctrl = None
        elif lower == "z":
            codes.append(MplPath.CLOSEPOLY)
            verts.append(tuple(start))
            cur = start.copy()
            last_cubic_ctrl = None
        else:
            raise ValueError(f"Unsupported SVG path command: {cmd}")

    return MplPath(verts, codes)


def draw_fontawesome_icon(ax, icon_name, x, y, size, color):
    icon = FA_ICONS[icon_name]
    view_w, view_h = icon["viewbox"]
    path = svg_path_to_mpl(icon["path"])
    scale = size / max(view_w, view_h)
    transform = (
        Affine2D()
        .scale(scale, -scale)
        .translate(x - view_w * scale / 2, y + view_h * scale / 2)
        + ax.transData
    )
    ax.add_patch(patches.PathPatch(path, facecolor=color, edgecolor="none", transform=transform, zorder=9))


def row_icon(ax, icon_name, row_x, row_y, row_w, row_h, color, badge_size=4.4, margin=1.4):
    x = row_x + row_w - margin - badge_size / 2
    y = row_y + row_h - margin - badge_size / 2
    rounded(
        ax,
        x - badge_size / 2,
        y - badge_size / 2,
        badge_size,
        badge_size,
        fc="white",
        ec=color,
        lw=1.0,
        r=0.75,
        z=8,
        alpha=0.92,
    )
    draw_fontawesome_icon(ax, icon_name, x, y, size=2.25, color=color)


def light_for(color):
    if color == PAL.green:
        return "#f5fbf7"
    if color == PAL.orange:
        return "#fff8f1"
    if color == PAL.azure:
        return "#f4fbff"
    if color == PAL.blue:
        return "#f6f9ff"
    return "#f8fafc"


def module(ax, x, y, w, h, title, ec=PAL.blue, fc="white", title_size=7.8):
    tint = light_for(ec)
    rounded(ax, x + 0.28, y - 0.28, w, h, fc="#142033", ec="none", lw=0, r=0.95, z=1, alpha=0.08)
    rounded(ax, x, y, w, h, fc=fc, ec=ec, lw=1.15, r=0.95, z=2)
    title_band = min(5.4, h * 0.36)
    rounded(
        ax,
        x + 0.25,
        y + h - title_band - 0.25,
        w - 0.5,
        title_band,
        fc=tint,
        ec="none",
        lw=0,
        r=0.65,
        z=3,
    )
    ax.plot(
        [x + 1.0, x + w - 1.0],
        [y + h - title_band - 0.35, y + h - title_band - 0.35],
        color=ec,
        lw=0.55,
        alpha=0.28,
        zorder=4,
    )
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
    content_x = x + 0.85
    content_y = y + 0.85
    content_w = w - 1.7
    content_h = h - title_band - 1.45
    rounded(
        ax,
        content_x,
        content_y,
        content_w,
        content_h,
        fc="#ffffff",
        ec="#eef2f7",
        lw=0.55,
        r=0.55,
        z=3,
        alpha=0.78,
    )
    return content_x + 0.35, content_y + 0.35, content_w - 0.7, content_h - 0.7


def group_card(ax, x, y, w, h, ec=PAL.blue):
    rounded(ax, x + 0.28, y - 0.28, w, h, fc="#142033", ec="none", lw=0, r=1.05, z=1, alpha=0.07)
    rounded(ax, x, y, w, h, fc="white", ec=ec, lw=1.15, r=1.05, z=2)
    rounded(
        ax,
        x + 0.35,
        y + 0.35,
        w - 0.7,
        h - 0.7,
        fc="#ffffff",
        ec="#eef2f7",
        lw=0.45,
        r=0.75,
        z=2,
        alpha=0.65,
    )


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


def draw_asset(ax, name, x, y, w, h, preserve_aspect=False):
    img = mpimg.imread(ASSETS[name])
    if preserve_aspect:
        img_h, img_w = img.shape[:2]
        img_aspect = img_w / img_h
        p0 = ax.transData.transform((0, 0))
        px = ax.transData.transform((1, 0))
        py = ax.transData.transform((0, 1))
        sx = abs(px[0] - p0[0])
        sy = abs(py[1] - p0[1])
        box_aspect = (w * sx) / (h * sy)
        if box_aspect > img_aspect:
            draw_h = h
            draw_w = (h * sy * img_aspect) / sx
        else:
            draw_w = w
            draw_h = (w * sx / img_aspect) / sy
        x = x + (w - draw_w) / 2
        y = y + (h - draw_h) / 2
        w, h = draw_w, draw_h
    ax.imshow(img, extent=(x, x + w, y, y + h), aspect="auto", zorder=4)


def build():
    canvas_h = 132.0
    fig, ax = plt.subplots(figsize=(12, 15.55), dpi=180)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, canvas_h)
    ax.axis("off")
    fig.patch.set_facecolor(PAL.bg)
    ax.set_facecolor(PAL.bg)

    row_x, row_w = 1.2, 97.6
    row_h = 26.0
    iter_row_h = 38.0
    card_h = 16.8
    row_gap = 4.0
    card_bottom = 2.3

    row1_y = 104.0
    row2_y = row1_y - row_h - row_gap
    row3_y = row2_y - iter_row_h - row_gap
    row4_y = row3_y - row_h - row_gap

    layer(ax, 1, "Inputs and dependencies", row_x, row1_y, row_w, row_h, PAL.azure, PAL.azure_light)
    layer(ax, 2, "Core PyDFT construction", row_x, row2_y, row_w, row_h, PAL.green, PAL.green_light)
    layer(ax, 3, "Kohn-Sham iteration", row_x, row3_y, row_w, iter_row_h, PAL.blue, PAL.blue_light)
    layer(ax, 4, "Exposed learning outputs", row_x, row4_y, row_w, row_h, PAL.orange, PAL.orange_light)
    row_icon(ax, "inbox", row_x, row1_y, row_w, row_h, PAL.azure)
    row_icon(ax, "gear", row_x, row2_y, row_w, row_h, PAL.green)
    row_icon(ax, "arrows-rotate", row_x, row3_y, row_w, iter_row_h, PAL.blue)
    row_icon(ax, "chart-line", row_x, row4_y, row_w, row_h, PAL.orange)

    # Layer 1
    left_x, left_y, left_w, left_h = 3.2, row1_y + card_bottom, 47.0, card_h
    right_x, right_y, right_w, right_h = 54.0, row1_y + card_bottom, 42.8, card_h
    group_card(ax, left_x, left_y, left_w, left_h, ec=PAL.azure)
    group_card(ax, right_x, right_y, right_w, right_h, ec=PAL.azure)
    ax.plot(
        [52.1, 52.1],
        [left_y + 2.4, left_y + left_h - 2.4],
        color=PAL.azure,
        lw=1.0,
        ls=(0, (4, 4)),
        alpha=0.32,
        zorder=3,
    )
    for i, (title, kind, sub) in enumerate(
        [
            ("Molecule", "molecule", "geometry"),
            ("Charge", "charge", "closed shell"),
            ("Basis set", "basis", "STO-3G, 6-31G(d)"),
        ]
    ):
        x = left_x + 2 + i * 15.0
        text(ax, x + 6.6, left_y + left_h - 2.3, title, size=10.4, weight="bold")
        draw_asset(ax, kind, x + 1.0, left_y + 4.35, 11.2, 8.4)
        text(ax, x + 6.6, left_y + 2.1, sub, size=7.0, color=PAL.ink)

    for i, (title, kind, sub) in enumerate(
        [
            ("PyQInt", "pyqint", "analytic integrals"),
            ("NumPy/SciPy", "numpy", "arrays, harmonics"),
            ("PyLebedev", "pylebedev", "angular quadrature"),
        ]
    ):
        x = right_x + 2.0 + i * 13.4
        text(ax, x + 5.8, right_y + right_h - 2.3, title, size=10.4, weight="bold", color=PAL.azure_dark)
        draw_asset(ax, kind, x + 0.7, right_y + 4.35, 10.2, 8.4)
        text(ax, x + 5.8, right_y + 2.2, sub, size=6.8, color=PAL.ink)

    arrow(ax, 50, row1_y - 0.3, 50, row2_y + row_h + 0.7, color=PAL.azure, lw=2.0)

    # Layer 2
    core_items = [
        ("Basis\namplitudes", "basis-amplitudes"),
        ("$S_{\\mu\\nu}, T_{\\mu\\nu}, V_{\\mu\\nu}$", "one-electron-integrals"),
        ("Becke grid\n$\\{\\mathbf{r}_i,w_i\\}$", "becke-grid"),
        ("Fuzzy weights\n$w_A(\\mathbf{r})$", "fuzzy-weights"),
        ("Lebedev\n$\\int d\\Omega$", "lebedev-grid"),
        ("Radial grid\n$r_k$", "radial-grid"),
    ]
    x0, y0, gap = 3.2, row2_y + card_bottom, 1.7
    mw = (93.6 - gap * 5) / 6
    for i, (title, kind) in enumerate(core_items):
        x = x0 + i * (mw + gap)
        cx, cy, cw, ch = module(ax, x, y0, mw, card_h, title, ec=PAL.green, title_size=8.45)
        draw_asset(ax, kind, cx, cy + 0.2, cw, ch - 0.2, preserve_aspect=True)

    arrow(ax, 50, row2_y - 0.1, 50, row3_y + iter_row_h + 0.3, color=PAL.green, lw=2.0)

    # Layer 3
    iter_items = [
        ("Density matrix\n$P_{\\mu\\nu}$", "density-matrix"),
        ("$\\rho(\\mathbf{r})$\n$\\nabla\\rho(\\mathbf{r})$", "density-gradient"),
        ("Hartree\n$v_H(\\mathbf{r})$", "hartree-potential"),
        ("XC\n$v_{xc}^{\\mathrm{LDA/PBE}}$", "xc-potential"),
        ("Build Fock\n$F=H+J+V_{xc}$", "fock-matrix-build"),
        ("Solve KS\n$FC=SC\\epsilon$", "solve-ks"),
        ("Update density\n$P=2CC^T$", "update-density"),
    ]
    ix0, iy0, igap = 3.2, row3_y + 11.5, 2.2
    iw = (93.6 - igap * 6) / 7
    centers = []
    cards = []
    for i, (title, kind) in enumerate(iter_items):
        x = ix0 + i * (iw + igap)
        cards.append((x, iy0, iw, card_h))
        cx, cy, cw, ch = module(ax, x, iy0, iw, card_h, title, ec=PAL.blue, title_size=7.9)
        draw_asset(ax, kind, cx, cy + 0.2, cw, ch - 0.2, preserve_aspect=True)
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
    diamond_y = row3_y + 6.6
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
        color="#c52020",
        lw=1.5,
        rad=9,
    )
    row4_top = row4_y + row_h
    arrow(ax, diamond_x, diamond_y - diamond_h / 2, diamond_x, row4_top + 0.1, color=PAL.green, lw=2.0)
    badge(ax, 35.0, diamond_y - 0.1, "No", "#c52020", fc="#fff5f5", w=4.8, h=2.0)
    yes_y = 0.5 * (diamond_y - diamond_h / 2 + row4_top + 0.1)
    badge(ax, diamond_x, yes_y, "Yes", PAL.green, fc="#f0faf4", w=5.2, h=2.0)

    # Layer 4
    output_items = [
        ("Matrices\n$F,P,S$", "output-matrices"),
        ("Energy terms\n$E_i$", "energy-terms"),
        ("Orbitals\n$\\psi_i$", "orbitals"),
        ("Density maps\n$\\rho(\\mathbf{r})$", "density-maps"),
        ("Becke cells\n$w_A(\\mathbf{r})$", "becke-cells"),
        ("Spherical coeff.\n$\\rho_{klm}$", "spherical-coefficients"),
        ("Timing\nbreakdown", "timing-breakdown"),
        ("Interactive\nexploration", "interactive-exploration"),
    ]
    ox0, oy0, ogap = 3.2, row4_y + card_bottom, 1.15
    ow = (93.6 - ogap * 7) / 8
    oh = card_h
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
            title_size=7.55,
        )
        draw_asset(ax, kind, cx, cy + 0.3, cw, ch - 0.3, preserve_aspect=True)

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
