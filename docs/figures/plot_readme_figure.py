#!/usr/bin/env python3
# Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology, Abu Dhabi, UAE)
"""Figure of the README: the head-on collision of two solitary waves computed by SerreGravityWave.m.

Run serre_readme_data.m first (GNU Octave), then, from this folder,
    python3 plot_readme_figure.py
which writes serre-collision.png. The style is the stamped copy manuscript_style.py beside it.
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import manuscript_style as hs  # noqa: E402

t_all, x_all, eta_all = np.loadtxt(HERE / "eta.csv", delimiter=",", unpack=True)
times = np.unique(t_all)
x = np.unique(x_all)
eta = eta_all.reshape(x.size, times.size).T          # eta[i, j]: time times[i], abscissa x[j]
crest = int(np.argmax(eta.max(axis=1)))               # the stored profile of largest elevation
shown = [0, crest, times.size - 1]
level = np.ceil(100 * np.abs(eta).max()) / 100

hs.configure(family="libertine")
plt.rcParams["savefig.dpi"] = 260
plt.rcParams["figure.facecolor"] = hs.PAPER
plt.rcParams["savefig.facecolor"] = hs.PAPER

fig = plt.figure(figsize=(hs.WIDTH, 4.0))
grid = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.025], height_ratios=[1.0, 1.5])
top = fig.add_subplot(grid[0, 0])
bottom = fig.add_subplot(grid[1, 0], sharex=top)
bar_axis = fig.add_subplot(grid[1, 1])

for n, i in enumerate(shown):
    style = hs.series(n)
    style["marker"] = None
    top.plot(x, eta[i], lw=1.1, label=rf"$t = {times[i]:g}$", **style)
top.axhline(0.0, color=hs.SLATE, lw=0.5, linestyle=(0, (1.6, 1.6)))
top.set_ylabel(r"$\eta/d$")
top.tick_params(labelbottom=False)
hs.panel(top, "(a) free surface before, during and after the collision")
hs.legend_outside(top, where="above", ncol=3)

mesh = bottom.pcolormesh(x, times, eta, shading="gouraud", cmap=hs.SIGNED_COLORMAP,
                         norm=Normalize(-level, level), rasterized=True)
for i in shown:
    bottom.axhline(times[i], color=hs.SLATE, lw=0.6, linestyle=(0, (1.6, 1.6)))
bottom.set_xlim(x.min(), x.max())
bottom.set_ylim(times.min(), times.max())
bottom.set_xlabel(r"$x/d$")
bottom.set_ylabel(r"$t\sqrt{g/d}$")
hs.panel(bottom, r"(b) elevation $\eta/d$ in the plane $(x, t)$; dotted lines: the instants of (a)")

bar = fig.colorbar(mesh, cax=bar_axis, label=r"$\eta/d$")
bar.outline.set_linewidth(0.5)
bar.outline.set_edgecolor(hs.INK)
bar.ax.tick_params(width=0.5, length=2.5)

out = HERE / "serre-collision.png"
hs.save(fig, out)
# The background is opaque: without its alpha channel, and compressed a second time without
# loss, the file is smaller.
Image.open(out).convert("RGB").save(out, optimize=True)
print("written", out, f"(colour scale from -{level:g} to {level:g}; largest elevation at t = {times[crest]:g})")
