# SPDX-License-Identifier: GPL-3.0-only
"""Shared plumbing for the atlas generators (matplotlib → SVG, text kept as text)."""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
plt.rcParams.update({"svg.fonttype": "none", "font.family": "monospace", "font.size": 8,
                     "axes.spines.top": False, "axes.spines.right": False})
INK, NORTH, SOUTH, HOT, GREY = "#1b1b1f", "#2a6fb0", "#c9732a", "#c0392b", "#8a8a96"


def save(fig, name):
    path = HERE / name
    fig.savefig(path, format="svg", bbox_inches="tight")
    plt.close(fig)
    return path


def view(xyz, elev=24.0, azim=-35.0):
    """Orthographic view of 3D points: returns (x, up, depth); depth > 0 is toward the viewer."""
    a, e = np.radians(azim), np.radians(elev)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(e), -np.sin(e)], [0, np.sin(e), np.cos(e)]])
    q = np.asarray(xyz, float) @ Rz.T @ Rx.T
    return q[:, 0], q[:, 2], -q[:, 1]
