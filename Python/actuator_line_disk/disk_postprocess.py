r"""
=======================================================================================
 ____          ____    __    ______     __________   __      __       __        __
 \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
      \    \  |    |   ________________________________________________________________
       \    \ |    |  |  ______________________________________________________________|
        \    \|    |  |  |         __          __     __     __     ______      _______
         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/

  This file is part of VirtualFluids. VirtualFluids is free software: you can
  redistribute it and/or modify it under the terms of the GNU General Public
  License as published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  SPDX-License-Identifier: GPL-3.0-or-later
  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

=======================================================================================
"""

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import fsolve
import vf_postproc.definitions as definitions
plt.style.use(["vf_postproc.vf_style"])

def load_velocities(path_planes, thrust_coefficients, velocity_inlet, disk_center, diameter):

    dict_vx = {"lbm": {}, "analytic_disk": {}, "analytic_far": {}}
    def _solve_a(ct):
        f = lambda a: 4.0 * a * (1.0 - a) - ct
        return float(fsolve(f, 0.1)[0])

    a = np.array([_solve_a(ct) for ct in thrust_coefficients], dtype=float)
    velocities_disk = (1.0 - a) * velocity_inlet
    velocities_far = (1.0 - 2.0 * a) * velocity_inlet

    for index, ct in enumerate(thrust_coefficients):
        path_planes_case = path_planes / f"ct_0p{int(ct*100)}" / "zplane"
        print(path_planes_case)
        plane_files = f"zPlane_bin_lev_*_ID_0_Part_1_t_*.vtk.bin.vtu"
        plane_file = max(
            path_planes_case.glob(plane_files),
            key=lambda p: (
                int(p.stem.split("_lev_")[1].split("_")[0]),
                int(p.stem.split("_t_")[1].split(".")[0])
            )
        )
        plane = pv.read(plane_file)

        pts = plane.points
        x = pts[:, 0]
        y = pts[:, 1]
        vx_mean = plane.point_data["vx_mean"]

        coord_y = np.round(y, 10)
        y_unique = np.unique(coord_y)
        y_line = y_unique[np.argmin(np.abs(y_unique - disk_center[1]))]
        mask_y = np.isclose(coord_y, y_line)

        # x_raw = (x[mask_y] - position_disk[0]) / disk_diameter
        x_raw = x[mask_y] / diameter
        vx_raw = vx_mean[mask_y] / velocity_inlet

        x_unique = np.unique(x_raw)
        vx_profile = np.array([vx_raw[x_raw == xx].mean() for xx in x_unique])

        idx = np.argsort(x_unique)
        x_profile = x_unique[idx]
        vx_profile = vx_profile[idx]

        key = f"{ct:.1f}"
        dict_vx["lbm"][key] = np.column_stack((x_profile, vx_profile))
        dict_vx["analytic_disk"][key] = velocities_disk[index] / velocity_inlet
        dict_vx["analytic_far"][key] = velocities_far[index] / velocity_inlet

    return dict_vx


def plot_velocities(dict_vx, path_out, disk_x):

    disk_x = disk_x[0][0]

    # ==========================================================
    # Plot settings
    # ==========================================================
    color_disk_position = "#000000"
    color_disk_velocity = "#332288"
    color_disk_far = "#332288"
    color_lbm = "#44AA99"

    linestyle_disk_position = "-"
    linestyle_disk_velocity = "-."
    linestyle_disk_far = "--"
    linestyle_lbm = "-"

    linewidth_disk_position = 0.8
    linewidth_disk_velocity = 0.8
    linewidth_disk_far = 0.8
    linewidth_lbm = 0.9

    label_disk = "Analytical (disk)"
    label_far = "Analytical (far wake)"
    label_lbm = "LBM (VirtualFluids)"

    xlim = [-10,10]
    ylim = [0.65, 1.05]
    xticks = [-10,0,10]

    aspect_velocities = 0.2
    aspect_subplot = 0.4

    xlegend, ylegend = 0.80, 0.50

    margin_left = 0.1
    margin_right = 0.85
    margin_bottom = 0.22
    margin_top = 0.82
    wspace = -0.2
    hspace = 0.0

    # ==========================================================
    # Figure / layout
    # ==========================================================
    ct_keys = list(dict_vx["lbm"].keys())

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(ct_keys),
        sharex=True,
        sharey=True,
        figsize=(definitions.FIGWIDTH, definitions.FIGWIDTH * aspect_velocities),
        facecolor="white"
    )

    if len(ct_keys) == 1:
        axes = [axes]

    fig.subplots_adjust(
        left=margin_left,
        right=margin_right,
        bottom=margin_bottom,
        top=margin_top,
        wspace=wspace,
        hspace=hspace
    )

    # ==========================================================
    # Panels
    # ==========================================================
    for i, (ax, ct_key) in enumerate(zip(axes, ct_keys)):
        arr = np.asarray(dict_vx["lbm"][ct_key], dtype=float)
        x = arr[:, 0]
        u = arr[:, 1]

        u_disk = float(dict_vx["analytic_disk"][ct_key])
        u_far = float(dict_vx["analytic_far"][ct_key])

        ax.axhline(
            u_disk,
            color=color_disk_velocity,
            linestyle=linestyle_disk_velocity,
            linewidth=linewidth_disk_velocity,
            label=label_disk if i == 0 else None
        )
        ax.axhline(
            u_far,
            color=color_disk_far,
            linestyle=linestyle_disk_far,
            linewidth=linewidth_disk_far,
            label=label_far if i == 0 else None
        )

        ax.plot(
            x - disk_x, u,
            color=color_lbm,
            linestyle=linestyle_lbm,
            linewidth=linewidth_lbm,
            label=label_lbm if i == 0 else None
        )
        ax.axvline(
            disk_x - disk_x,
            color=color_disk_position,
            linestyle=linestyle_disk_position,
            linewidth=linewidth_disk_position
        )

        ax.set_title(rf"$C_T = {ct_key}$")
        ax.text(disk_x-disk_x/2, 1.045, "disk position", ha="center", va="top", color=color_disk_position)
        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_ylim(ylim)
        ax.set_xlabel(r"$x/D$")
        ax.set_box_aspect(aspect_subplot)

        if i > 0:
            ax.tick_params(axis="y", left=False, labelleft=False)

    axes[0].set_ylabel(r"$u/u_\infty$")

    # ==========================================================
    # Legend
    # ==========================================================
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", ncol=1, frameon=False, bbox_to_anchor=(xlegend, ylegend))

    plt.savefig(path_out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)