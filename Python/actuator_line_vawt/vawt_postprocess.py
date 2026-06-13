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
#!/usr/bin/env python3

import numpy as np
import pyvista as pv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vf_postproc.definitions as definitions
import os, re
from collections import defaultdict
plt.style.use(["vf_postproc.vf_style"])

# =========================
# Parameters
# =========================
x_over_D_list = [0.75,1,1.25,1.5,1.75,2,5,10]
x_over_D_validation = x_over_D_list[:6]

Uinf = 9.3
z_target = 4
y_target = 4
D_rotor = 1.0
R = D_rotor / 2.0

tol = 1.5
smooth = True
n_interp = 200

PATH_SCRIPT = Path(__file__).resolve().parent
PATH_PLANES = PATH_SCRIPT / "data" / "output" / "planes"
PATH_REFERENCE = PATH_SCRIPT / "data" / "validation" / "wake"
PATH_POST = PATH_SCRIPT / "data" / "post"
if not os.path.exists(PATH_POST): os.makedirs(PATH_POST)

eps_list_lbm = []
eps_list_fvm = []
vx_y_dict = {"exp":{},"fvm":{},"lbm":{},"lbm_std":{}}
vx_z_dict = {"exp":{},"fvm":{},"lbm":{},"lbm_std":{}}
vy_y_dict = {"exp":{},"fvm":{},"lbm":{},"lbm_std":{}}

color_experiment = "#332288"
color_fvm = "#CC6677"
color_lbm = "#44AA99"

linestyle_exp = "None"
linestyle_fvm = "--"
linestyle_lbm = "-"

marker_exp = "+"
marker_fvm = "None"
marker_lbm = "None"

linewidth_exp = 0.1
linewidth_fvm = 0.9
linewidth_lbm = linewidth_fvm

markersize_exp = 6.0
markersize_fvm = 0.0
markersize_lbm = 0.0

alpha_exp = 1.0
alpha_fvm = 1.0
alpha_lbm = 1.0
color_lbm_std_alpha = 0.4

def load_wake(dx_fine):

    for index, x_over_D in enumerate(x_over_D_list):

        # ==========================================================
        # SIMULATION
        # ==========================================================
        plane_files = f"x_plane_{index}_bin_lev_*_ID_0_Part_1_t_*.vtk.bin.vtu"
        plane_file = max(PATH_PLANES.glob(plane_files),key=lambda p: (
                        int(p.stem.split("_lev_")[1].split("_")[0]), # finest level
                        int(p.stem.split("_t_")[1].split(".")[0]))) # last timestep
        plane = pv.read(PATH_PLANES / plane_file)

        pts = plane.points
        y = pts[:, 1]
        z = pts[:, 2]
        vx_mean = plane.point_data["vx_mean"]
        vx_std = np.sqrt(plane.point_data["vx_var"])
        vy_mean = plane.point_data["vy_mean"]
        vy_std = np.sqrt(plane.point_data["vy_var"])

        # ==========================================================
        # LBM (VirtualFluids)
        # ==========================================================

        coord_z = np.round(z, 10)
        z_unique = np.unique(coord_z)
        z_line = z_unique[np.argmin(np.abs(z_unique - z_target))]
        mask_z = np.isclose(coord_z, z_line)

        # Horizontal ux profile
        vx_y_raw_mean = vx_mean[mask_z] / Uinf
        vx_y_raw_std = vx_std[mask_z] / Uinf
        y_prof = (y[mask_z] - y_target) / D_rotor

        y_unique = np.unique(y_prof)
        vx_y_mean = np.array([vx_y_raw_mean[y_prof == yy].mean() for yy in y_unique])
        vx_y_std = np.array([vx_y_raw_std[y_prof == yy].mean() for yy in y_unique])

        idx = np.argsort(y_unique)
        y_sim = y_unique[idx]
        vx_y_mean = vx_y_mean[idx]
        vx_y_std = vx_y_std[idx]

        # Vertical ux profile

        coord_y = np.round(y, 10)
        y_unique_plane = np.unique(coord_y)
        y_line = y_unique_plane[np.argmin(np.abs(y_unique_plane - y_target))]
        mask_y = np.isclose(coord_y, y_line)

        vx_z_raw_mean = vx_mean[mask_y] / Uinf
        vx_z_raw_std = vx_std[mask_y] / Uinf
        z_prof = (z[mask_y] - z_target) / D_rotor

        z_unique = np.unique(z_prof)
        vx_z_mean = np.array([vx_z_raw_mean[z_prof == zz].mean() for zz in z_unique])
        vx_z_std = np.array([vx_z_raw_std[z_prof == zz].mean() for zz in z_unique])

        idx = np.argsort(z_unique)
        z_sim = z_unique[idx]
        vx_z_mean = vx_z_mean[idx]
        vx_z_std = vx_z_std[idx]

        # Horizontal uy profile

        vy_y_raw_mean = vy_mean[mask_z] / Uinf
        vy_y_raw_std = vy_std[mask_z] / Uinf

        vy_y_mean = np.array([vy_y_raw_mean[y_prof == yy].mean() for yy in y_unique])
        vy_y_std = np.array([vy_y_raw_std[y_prof == yy].mean() for yy in y_unique])

        idx = np.argsort(y_unique)
        vy_y_mean = vy_y_mean[idx]
        vy_y_std = vy_y_std[idx]

        # collect LBM data in Nx2 shape

        vx_y_dict["lbm"][f"{x_over_D*100}"] = np.column_stack((y_sim,vx_y_mean))
        vx_z_dict["lbm"][f"{x_over_D*100}"] = np.column_stack((z_sim,vx_z_mean))
        vy_y_dict["lbm"][f"{x_over_D*100}"] = np.column_stack((y_sim,vy_y_mean))
        vx_y_dict["lbm_std"][f"{x_over_D*100}"] = np.column_stack((y_sim,vx_y_std))
        vx_z_dict["lbm_std"][f"{x_over_D*100}"] = np.column_stack((z_sim,vx_z_std))
        vy_y_dict["lbm_std"][f"{x_over_D*100}"] = np.column_stack((y_sim,vy_y_std))

        # ==========================================================
        # FVM (OpenFoam)
        # ==========================================================

        vx_y_fvm = np.loadtxt(PATH_REFERENCE / "Rogowski_2025" / f"ux_y_{int(x_over_D*100)}o100.dat",comments="%")
        vx_z_fvm = np.loadtxt(PATH_REFERENCE / "Rogowski_2025" / f"ux_z_{int(x_over_D*100)}o100.dat",comments="%")
        vy_y_fvm = np.loadtxt(PATH_REFERENCE / "Rogowski_2025" / f"uy_y_{int(x_over_D*100)}o100.dat",comments="%")

        vx_y_dict["fvm"][f"{x_over_D*100}"] = vx_y_fvm
        vx_z_dict["fvm"][f"{x_over_D*100}"] = vx_z_fvm
        vy_y_dict["fvm"][f"{x_over_D*100}"] = vy_y_fvm

        # ==========================================================
        # Experimental (PIV)
        # ==========================================================

        # if index < 6: # left out because of missing license
        #     vx_y_exp = np.loadtxt(PATH_REFERENCE / "Tescione_2014" / f"ux_y_{int(x_over_D*100)}o100.dat")
        #     vx_z_exp = np.loadtxt(PATH_REFERENCE / "Mendoza_2019" / f"ux_z_{int(x_over_D*100)}o100.dat")
        #     vy_y_exp = np.loadtxt(PATH_REFERENCE / "Tescione_2014" / f"uy_y_{int(x_over_D*100)}o100.dat")

        #     vx_y_dict["exp"][f"{x_over_D*100}"] = vx_y_exp
        #     vx_z_dict["exp"][f"{x_over_D*100}"] = vx_z_exp
        #     vy_y_dict["exp"][f"{x_over_D*100}"] = vy_y_exp

    return vx_y_dict, vx_z_dict, vy_y_dict

def compute_wake_error(vx_y_dict, vx_z_dict, vy_y_dict):

    def _calculate_error(v_ref,v):
        dv_ref = v_ref[:] - 1
        dv = v[:] - 1
        return np.sqrt(np.sum((dv_ref - dv)**2) / np.sum(dv_ref**2))
    
    for index, x_over_D in enumerate(x_over_D_validation):

        vx_y_exp = vx_y_dict["exp"][f"{x_over_D*100}"]
        vx_z_exp = vx_z_dict["exp"][f"{x_over_D*100}"]
        vy_y_exp = vy_y_dict["exp"][f"{x_over_D*100}"]

        vx_y_fvm = vx_y_dict["fvm"][f"{x_over_D*100}"]
        vx_z_fvm = vx_z_dict["fvm"][f"{x_over_D*100}"]
        vy_y_fvm = vy_y_dict["fvm"][f"{x_over_D*100}"]

        vx_y_lbm = vx_y_dict["lbm"][f"{x_over_D*100}"]
        vx_z_lbm = vx_z_dict["lbm"][f"{x_over_D*100}"]
        vy_y_lbm = vy_y_dict["lbm"][f"{x_over_D*100}"]

        vx_y_lbm_std = vx_y_dict["lbm_std"][f"{x_over_D*100}"]
        vx_z_lbm_std = vx_z_dict["lbm_std"][f"{x_over_D*100}"]
        vy_y_lbm_std = vy_y_dict["lbm_std"][f"{x_over_D*100}"]
    
        # ==========================================================
        # error
        # ==========================================================

        vx_y_lbm_interp = np.interp(vx_y_exp[:,0], vx_y_lbm[:,0], vx_y_lbm[:,1])
        vx_y_fvm_interp = np.interp(vx_y_exp[:,0], vx_y_fvm[:,0], vx_y_fvm[:,1])

        z_min = max(vx_z_exp[:,0].min(), vx_z_lbm[:,0].min(), vx_z_fvm[:,0].min())
        z_max = min(vx_z_exp[:,0].max(), vx_z_lbm[:,0].max(), vx_z_fvm[:,0].max())
        mask_z = (vx_z_exp[:,0] >= z_min) & (vx_z_exp[:,0] <= z_max)
        vx_z_lbm_interp = np.interp(vx_z_exp[mask_z,0], vx_z_lbm[:,0], vx_z_lbm[:,1])
        vx_z_fvm_interp = np.interp(vx_z_exp[mask_z,0], vx_z_fvm[:,0], vx_z_fvm[:,1])

        err_vx_y_lbm = _calculate_error(vx_y_exp[:,1],vx_y_lbm_interp)
        err_vx_y_fvm = _calculate_error(vx_y_exp[:,1],vx_y_fvm_interp)

        err_vx_z_lbm = _calculate_error(vx_z_exp[mask_z,1],vx_z_lbm_interp)
        err_vx_z_fvm = _calculate_error(vx_z_exp[mask_z,1],vx_z_fvm_interp)

        # ==========================================================
        # final error
        # ==========================================================

        epsilon_lbm = 0.5 * (err_vx_y_lbm + err_vx_z_lbm)
        epsilon_fvm = 0.5 * (err_vx_y_fvm + err_vx_z_fvm)
        eps_list_lbm.append(epsilon_lbm)
        eps_list_fvm.append(epsilon_fvm)

    epsilon_total_lbm = np.sum(eps_list_lbm)
    epsilon_total_fvm = np.sum(eps_list_fvm)

    return epsilon_total_lbm, epsilon_total_fvm

def plot_wake(vx_y_dict, vx_z_dict, vy_y_dict, path_out):

    import matplotlib.pyplot as plt
    import numpy as np
    from pathlib import Path
    PATH_SCRIPT = Path(__file__).resolve().parent
    PATH_POST = PATH_SCRIPT / "data" / "post"

    # ==========================================================
    # Plot settings
    # ==========================================================

    label_exp = "Experiment (PIV)"
    label_fvm = "FVM (OpenFOAM)"
    label_lbm = "LBM (VirtualFluids)"

    xlim_vx = [0.15,1.05]
    xlim_vy = [-0.1,0.15]
    ylim_y = [-0.9,0.9]
    ylim_z = [0,0.55]

    aspect_wake = 0.6

    xlegend,ylegend = 0.8,0.5

    # margins / spacing
    margin_left = 0.12
    margin_right = 0.8
    margin_bottom = 0.02
    margin_top = 0.95
    wspace = 0.0
    hspace = 0.0
    gap_height_between_row_2_3 = 0.6

    # ==========================================================
    # Data / layout
    # ==========================================================
    x_keys = list(vx_y_dict["lbm"].keys())
    ncols = len(x_keys)
    import vf_postproc.definitions as definitions
    fig, axes = plt.subplots(
        nrows=4,
        ncols=ncols,
        sharex=False,
        sharey=False,
        figsize=(definitions.FIGWIDTH, definitions.FIGWIDTH * aspect_wake),
        gridspec_kw={"height_ratios": [1, 1, gap_height_between_row_2_3, 1]}
    )

    fig.subplots_adjust(
        left=margin_left,
        right=margin_right,
        bottom=margin_bottom,
        top=margin_top,
        wspace=wspace,
        hspace=hspace
    )

    # make sure axes is always 2D
    axes = np.atleast_2d(axes)
    for j in range(ncols):
        axes[0, j].tick_params(labelbottom=False)
        axes[1, j].sharex(axes[0, j])
        # hide x ticks on row 1 (shared)
        axes[0, j].tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        # hide y ticks for all but first column
        if j > 0:
            axes[0, j].tick_params(axis="y", which="both", left=False, labelleft=False)
            axes[1, j].tick_params(axis="y", which="both", left=False, labelleft=False)
            axes[3, j].tick_params(axis="y", which="both", left=False, labelleft=False)
        # define who shares and axis
        axes[0, j].sharey(axes[0, 0])   # row 1 shares y/D
        axes[3, j].sharey(axes[3, 0])   # row 4 shares y/D
        axes[1, j].sharey(axes[1, 0])   # row 2 shares z/D only within row 2

    # ==========================================================
    # Loop over downstream planes
    # ==========================================================
    for j, key in enumerate(x_keys):

        # ------------------------------------------------------
        # Row 1: vx_y
        # x-axis = normalized velocity, y-axis = coordinate
        # ------------------------------------------------------
        ax = axes[0, j]

        # if key in vx_y_dict["exp"]:
        #     arr = np.asarray(vx_y_dict["exp"][key], dtype=float)
        #     ax.scatter(
        #         arr[:, 1], arr[:, 0],
        #         color=color_experiment,
        #         marker=marker_exp,
        #         s=markersize_exp,
        #         alpha=alpha_exp,
        #         label=label_exp if j == 0 else None
        #     )

        if key in vx_y_dict["fvm"]:
            arr = np.asarray(vx_y_dict["fvm"][key], dtype=float)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_fvm,
                linestyle=linestyle_fvm,
                marker=marker_fvm,
                linewidth=linewidth_fvm,
                markersize=markersize_fvm,
                alpha=alpha_fvm,
                label=label_fvm if j == 0 else None
            )

        if key in vx_y_dict["lbm"]:
            arr = np.asarray(vx_y_dict["lbm"][key], dtype=float)
            arr_std = np.asarray(vx_y_dict["lbm_std"][key], dtype=float)
            ax.fill_betweenx(arr[:,0],arr[:,1] - arr_std[:,1],arr[:,1] + arr_std[:,1],color=color_lbm,alpha=color_lbm_std_alpha,linewidth=0.0)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_lbm,
                linestyle=linestyle_lbm,
                marker=marker_lbm,
                linewidth=linewidth_lbm,
                markersize=markersize_lbm,
                alpha=alpha_lbm,
                label=label_lbm if j == 0 else None
            )

        if j == 0: 
            ax.set_title(rf"$x/D = {float(key)/100:.2f}$", loc='right')
        else : 
            ax.set_title(rf"${float(key)/100:.2f}$")

        ax.set_xlim(xlim_vx)
        ax.set_ylim(ylim_y)

        # ------------------------------------------------------
        # Row 2: vx_z
        # ------------------------------------------------------
        ax = axes[1, j]

        # if key in vx_z_dict["exp"]:
        #     arr = np.asarray(vx_z_dict["exp"][key], dtype=float)
        #     ax.scatter(
        #         arr[:, 1], arr[:, 0],
        #         color=color_experiment,
        #         marker=marker_exp,
        #         s=markersize_exp,
        #         alpha=alpha_exp,
        #         label=label_exp if j == 0 else None
        #     )

        if key in vx_z_dict["fvm"]:
            arr = np.asarray(vx_z_dict["fvm"][key], dtype=float)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_fvm,
                linestyle=linestyle_fvm,
                marker=marker_fvm,
                linewidth=linewidth_fvm,
                markersize=markersize_fvm,
                alpha=alpha_fvm
            )

        if key in vx_z_dict["lbm"]:
            arr = np.asarray(vx_z_dict["lbm"][key], dtype=float)
            arr_std = np.asarray(vx_z_dict["lbm_std"][key], dtype=float)
            ax.fill_betweenx(arr[:,0],arr[:,1] - arr_std[:,1],arr[:,1] + arr_std[:,1],color=color_lbm,alpha=color_lbm_std_alpha,linewidth=0.0)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_lbm,
                linestyle=linestyle_lbm,
                marker=marker_lbm,
                linewidth=linewidth_lbm,
                markersize=markersize_lbm,
                alpha=alpha_lbm
            )

        ax.set_xlim(xlim_vx)
        ax.set_ylim(ylim_z)

        # ------------------------------------------------------
        # Row 2/3 (just for x axis label)
        # ------------------------------------------------------
        ax = axes[2, j]
        ax.axis("off")

        # ------------------------------------------------------
        # Row 3: vy_y
        # ------------------------------------------------------
        ax = axes[3, j]

        # if key in vy_y_dict["exp"]:
        #     arr = np.asarray(vy_y_dict["exp"][key], dtype=float)
        #     ax.scatter(
        #         arr[:, 1], arr[:, 0],
        #         color=color_experiment,
        #         marker=marker_exp,
        #         s=markersize_exp,
        #         alpha=alpha_exp,
        #         label=label_exp if j == 0 else None
        #     )

        if key in vy_y_dict["fvm"]:
            arr = np.asarray(vy_y_dict["fvm"][key], dtype=float)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_fvm,
                linestyle=linestyle_fvm,
                marker=marker_fvm,
                linewidth=linewidth_fvm,
                markersize=markersize_fvm,
                alpha=alpha_fvm
            )

        if key in vy_y_dict["lbm"]:
            arr = np.asarray(vy_y_dict["lbm"][key], dtype=float)
            arr_std = np.asarray(vy_y_dict["lbm_std"][key], dtype=float)
            ax.fill_betweenx(arr[:,0],arr[:,1] - arr_std[:,1],arr[:,1] + arr_std[:,1],color=color_lbm,alpha=color_lbm_std_alpha,linewidth=0.0)
            ax.plot(
                arr[:, 1], arr[:, 0],
                color=color_lbm,
                linestyle=linestyle_lbm,
                marker=marker_lbm,
                linewidth=linewidth_lbm,
                markersize=markersize_lbm,
                alpha=alpha_lbm
            )

        ax.set_xlim(xlim_vy)
        ax.set_ylim(ylim_y)
        ax.tick_params(axis="x",labelrotation=45)
        ax.set_xticks([0,0.1])

    # ==========================================================
    # Labels
    # ==========================================================
    axes[0, 0].set_ylabel(r"$y/D$")
    axes[1, 0].set_ylabel(r"$z/D$")
    axes[3, 0].set_ylabel(r"$y/D$")

    for j in range(ncols):
        axes[1, j].set_xlabel(r"$u_x/u_\infty$")
        axes[3, j].set_xlabel(r"$u_y/u_\infty$")

    # ==========================================================
    # Legend
    # ==========================================================
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", ncol=1, frameon=False, bbox_to_anchor=(xlegend,ylegend))

    plt.savefig(path_out, dpi=300, facecolor = 'white')

    return fig, axes


def load_loads(number_blades, flag_flow_curvature, flag_end_effects):

    alpha_dict = {"fvm":{}, "lbm":{}}
    fn_dict = {"fvm":{}, "lbm":{}}
    ft_dict = {"fvm":{}, "lbm":{}}

    LOAD_FILE_PATTERN = re.compile(r"^ALM_ID_(?P<process_id>\d+)_t_(?P<t>\d+)\.bin\.vtu$")

    def parse_load_info(path_load):
        groups = LOAD_FILE_PATTERN.match(path_load.name).groupdict()
        return {"process_id": int(groups["process_id"]), "t": int(groups["t"])}

    def _sort_by_azimuth(x, y):
        x = np.mod(np.asarray(x, dtype=float), 360.0)
        order = np.argsort(x)
        return np.column_stack((x[order], np.asarray(y, dtype=float)[order]))
    
    def _azimuthal_mean_std(theta_deg, values, nbins=200):
        theta = np.mod(np.asarray(theta_deg, dtype=float), 360.0)
        values = np.asarray(values, dtype=float)
        edges = np.linspace(0.0, 360.0, nbins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        idx = np.digitize(theta, edges) - 1
        idx[idx == nbins] = 0

        sums = np.zeros(nbins)
        sums_sq = np.zeros(nbins)
        counts = np.zeros(nbins)
        np.add.at(sums, idx, values)
        np.add.at(sums_sq, idx, values * values)
        np.add.at(counts, idx, 1.0)

        mean = np.divide(sums, counts, out=np.full(nbins, np.nan), where=counts > 0.0)
        var = np.divide(sums_sq, counts, out=np.full(nbins, np.nan), where=counts > 0.0) - mean * mean
        std = np.sqrt(np.clip(var, 0.0, None))
        valid = counts > 0.0
        return centers[valid], np.column_stack((mean[valid], std[valid]))

    validation_dir = PATH_SCRIPT / "data" / "validation" / "loads"
    loads_dir = PATH_SCRIPT / "data" / "output" / "loads"

    # ==========================================================
    # FVM / reference
    # ==========================================================
    if (not flag_flow_curvature) and (not flag_end_effects):
        alpha_raw = np.loadtxt(validation_dir / "EllipSys3D" / "EllipSys3D_aoa.dat",comments="%")
        fn_raw = np.loadtxt(validation_dir / "EllipSys3D" / "EllipSys3D_fn.dat",comments="%")
        ft_raw = np.loadtxt(validation_dir / "EllipSys3D" / "EllipSys3D_ft.dat",comments="%")

        theta_alpha = np.mod(alpha_raw[:, 1], 360.0)
        theta_fn = np.mod(fn_raw[:, 1], 360.0)
        theta_ft = np.mod(ft_raw[:, 1], 360.0)

        ref_alpha = alpha_raw[:, alpha_raw.shape[1] // 2]
        ref_fn = fn_raw[:, fn_raw.shape[1] // 2]
        ref_ft = ft_raw[:, ft_raw.shape[1] // 2]
    else:
        alpha_raw = np.loadtxt(validation_dir / "Rogowski_2025" / "Rogowski_2025_alpha.dat",comments="%")
        fn_raw = np.loadtxt(validation_dir / "Rogowski_2025" / "Rogowski_2025_fn.dat",comments="%")
        ft_raw = np.loadtxt(validation_dir / "Rogowski_2025" / "Rogowski_2025_ft.dat",comments="%")

        theta_alpha = np.mod(alpha_raw[:, 0], 360.0)
        theta_fn = np.mod(fn_raw[:, 0], 360.0)
        theta_ft = np.mod(ft_raw[:, 0], 360.0)

        ref_alpha = alpha_raw[:, 1]
        ref_fn = fn_raw[:, 1]
        ref_ft = ft_raw[:, 1]

    alpha_dict["fvm"] = _sort_by_azimuth(theta_alpha, ref_alpha)
    fn_dict["fvm"] = _sort_by_azimuth(theta_fn, ref_fn)
    ft_dict["fvm"] = _sort_by_azimuth(theta_ft, ref_ft)

    # ==========================================================
    # LBM (VirtualFluids)
    # ==========================================================
    grouped = defaultdict(list)
    for fp in sorted(loads_dir.glob("ALM_ID_*_t_*.bin.vtu")):
        info = parse_load_info(fp)
        grouped[info["t"]].append((info["process_id"], fp))

    theta_vals, alpha_vals, fn_vals, ft_vals = [], [], [], []
    for timestep in sorted(grouped):
        arrays = [pv.read(fp) for _, fp in sorted(grouped[timestep], key=lambda item: item[0])]
        theta = np.concatenate([arr.point_data["azimuthDeg"].ravel() for arr in arrays])
        alpha = np.concatenate([arr.point_data["angleOfAttackDeg"].ravel() for arr in arrays])
        fn = np.concatenate([arr.point_data["forceNormal"].ravel() for arr in arrays])
        ft = np.concatenate([arr.point_data["forceTangential"].ravel() for arr in arrays])

        center_index = (theta.size // number_blades) // 2
        theta_vals.append(np.mod(theta[center_index], 360.0))
        alpha_vals.append(alpha[center_index])
        fn_vals.append(fn[center_index])
        ft_vals.append(ft[center_index])

    lbm_theta, lbm_alpha = _azimuthal_mean_std(theta_vals, alpha_vals)
    _, lbm_fn = _azimuthal_mean_std(theta_vals, fn_vals)
    _, lbm_ft = _azimuthal_mean_std(theta_vals, ft_vals)

    alpha_dict["lbm"] = _sort_by_azimuth(lbm_theta, lbm_alpha)
    fn_dict["lbm"] = _sort_by_azimuth(lbm_theta, lbm_fn)
    ft_dict["lbm"] = _sort_by_azimuth(lbm_theta, lbm_ft)

    return alpha_dict, fn_dict, ft_dict

def plot_loads(alpha_dict, fn_dict, ft_dict, flag_flow_curvature, flag_end_effects, path_out):

    # def plot_loads(alpha_dict, fn_dict, ft_dict, path_out):

    # path_out = Path(path_out)
    # path_out.parent.mkdir(parents=True, exist_ok=True)

    # ==========================================================
    # Plot settings
    # ==========================================================

    if flag_flow_curvature == 0 and flag_end_effects == 0 :
        label_fvm = "FVM (EllipSys3D)"
    else:
        label_fvm = "FVM (OpenFoam)"
    label_lbm = "LBM (VirtualFluids)"

    xlim_theta = [0.0, 360.0]
    xticks_theta = [0, 90, 180, 270, 360]

    aspect_loads = 0.4
    aspect_subplot_loads = 0.8

    xlegend, ylegend = 0.8, 0.4

    # margins / spacing
    margin_left = 0.02
    margin_right = 0.8
    margin_bottom = 0.18
    margin_top = 0.7
    wspace = 0.7
    hspace = 1

    # ==========================================================
    # Figure / layout
    # ==========================================================
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        sharex=True,
        sharey=False,
        figsize=(definitions.FIGWIDTH, definitions.FIGWIDTH * aspect_loads),
        facecolor="white"
    )

    fig.subplots_adjust(
        left=margin_left,
        right=margin_right,
        bottom=margin_bottom,
        top=margin_top,
        wspace=wspace,
        hspace=hspace
    )

    # ==========================================================
    # Panel 1: alpha
    # ==========================================================
    ax = axes[0]

    arr = np.asarray(alpha_dict["fvm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_fvm,
        linestyle=linestyle_fvm,
        marker=marker_fvm,
        linewidth=linewidth_fvm,
        markersize=markersize_fvm,
        alpha=alpha_fvm,
        label=label_fvm
    )

    arr = np.asarray(alpha_dict["lbm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_lbm,
        linestyle=linestyle_lbm,
        marker=marker_lbm,
        linewidth=linewidth_lbm,
        markersize=markersize_lbm,
        alpha=alpha_lbm,
        label=label_lbm
    )

    ax.set_ylabel(r"$\alpha$ in deg")
    ax.set_xlim(xlim_theta)
    ax.set_xticks(xticks_theta)
    ax.set_xlabel(r"$\theta$ in deg")
    ax.set_box_aspect(aspect_subplot_loads)

    # ==========================================================
    # Panel 2: fn
    # ==========================================================
    ax = axes[1]

    arr = np.asarray(fn_dict["fvm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_fvm,
        linestyle=linestyle_fvm,
        marker=marker_fvm,
        linewidth=linewidth_fvm,
        markersize=markersize_fvm,
        alpha=alpha_fvm
    )

    arr = np.asarray(fn_dict["lbm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_lbm,
        linestyle=linestyle_lbm,
        marker=marker_lbm,
        linewidth=linewidth_lbm,
        markersize=markersize_lbm,
        alpha=alpha_lbm
    )

    ax.set_ylabel(r"$F_n/(q_0 R)$")
    ax.set_xlim(xlim_theta)
    ax.set_xticks(xticks_theta)
    ax.set_xlabel(r"$\theta$ in deg")
    ax.set_box_aspect(aspect_subplot_loads)

    # ==========================================================
    # Panel 3: ft
    # ==========================================================
    ax = axes[2]

    arr = np.asarray(ft_dict["fvm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_fvm,
        linestyle=linestyle_fvm,
        marker=marker_fvm,
        linewidth=linewidth_fvm,
        markersize=markersize_fvm,
        alpha=alpha_fvm
    )

    arr = np.asarray(ft_dict["lbm"], dtype=float)
    ax.plot(
        arr[:, 0], arr[:, 1],
        color=color_lbm,
        linestyle=linestyle_lbm,
        marker=marker_lbm,
        linewidth=linewidth_lbm,
        markersize=markersize_lbm,
        alpha=alpha_lbm
    )

    ax.set_ylabel(r"$F_t/(q_0 R)$")
    ax.set_xlim(xlim_theta)
    ax.set_xticks(xticks_theta)
    ax.set_xlabel(r"$\theta$ in deg")
    ax.set_box_aspect(aspect_subplot_loads)

    # ==========================================================
    # Legend
    # ==========================================================
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", ncol=1, frameon=False, bbox_to_anchor=(xlegend, ylegend))

    plt.savefig(path_out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)

    return alpha_dict, fn_dict, ft_dict


import subprocess, time, shutil, threading
class CPUEnergyMonitor:
    def __init__(self, path="/sys/class/powercap/intel-rapl:0/energy_uj", allow_sudo_fallback=True):
        self.path = Path(path)
        self.allow_sudo_fallback = allow_sudo_fallback

    def _read_energy_uj(self):
        try:
            return float(self.path.read_text().strip())
        except PermissionError:
            if not self.allow_sudo_fallback:
                raise
            return float(
                subprocess.check_output(
                    [shutil.which("sudo"), "-n", "cat", str(self.path)],
                    text=True,
                ).strip()
            )

    def start(self):
        self.t0 = time.perf_counter()
        self.e0 = self._read_energy_uj()

    def stop(self):
        dt = time.perf_counter() - self.t0
        de = self._read_energy_uj() - self.e0
        return {
            "duration_s": dt,
            "energy_wh": de / 3_600_000_000.0,
            "avg_power_w": (de / 1_000_000.0) / dt,
        }


class GPUEnergyMonitor:
    def __init__(self, gpu_index=0, poll_interval=1.0):
        self.gpu_index = str(gpu_index)
        self.poll_interval = poll_interval

    def _loop(self):
        while self.running:
            self.samples.append(
                float(
                    subprocess.check_output(
                        ["nvidia-smi", "--query-gpu=power.draw", "--format=csv,noheader,nounits", "-i", self.gpu_index],
                        text=True,
                    ).strip().splitlines()[0]
                )
            )
            time.sleep(self.poll_interval)

    def start(self):
        self.samples = []
        self.t0 = time.perf_counter()
        self.running = True
        self.thread = threading.Thread(target=self._loop, daemon=True)
        self.thread.start()

    def stop(self):
        self.running = False
        self.thread.join()
        dt = time.perf_counter() - self.t0
        avg_power_w = np.mean(self.samples)
        return {
            "duration_s": dt,
            "energy_wh": avg_power_w * dt / 3600.0,
            "avg_power_w": avg_power_w,
        }
    
def plot_performance(path_out, performance_fvm, performance_lbm):

    path_out = Path(path_out)
    path_out.parent.mkdir(parents=True, exist_ok=True)

    def as_two_values(values, name):
        if len(values) != 2:
            raise ValueError(f"'{name}' must contain exactly two values: [time_per_rev, energy_per_rev].")
        return values[0], values[1]

    def as_positive_float(value, name):
        value = float(value)
        if value <= 0.0:
            raise ValueError(f"'{name}' must be > 0 for a log-scale plot, got {value}.")
        return value

    fvm_time_raw, fvm_energy_raw = as_two_values(performance_fvm, "performance_fvm")
    lbm_time_raw, lbm_energy_raw = as_two_values(performance_lbm, "performance_lbm")

    time_fvm = as_positive_float(fvm_time_raw, "performance_fvm[0]")
    energy_fvm = as_positive_float(fvm_energy_raw, "performance_fvm[1]")
    time_lbm = as_positive_float(lbm_time_raw, "performance_lbm[0]")

    if np.isscalar(lbm_energy_raw):
        energy_lbm_gpu = as_positive_float(lbm_energy_raw, "performance_lbm[1]")
        energy_lbm_cpu = 0.0
    else:
        lbm_parts = np.asarray(list(lbm_energy_raw), dtype=float).ravel()
        if lbm_parts.size == 0 or lbm_parts.size > 2:
            raise ValueError("'performance_lbm[1]' must be a scalar or [gpu_energy, cpu_energy].")
        energy_lbm_gpu = as_positive_float(lbm_parts[0], "performance_lbm[1][0]")
        energy_lbm_cpu = as_positive_float(lbm_parts[1], "performance_lbm[1][1]") if lbm_parts.size == 2 else 0.0

    labels = ["FVM", "LBM"]
    x = np.arange(len(labels))
    width = 0.60
    aspect_performance = 0.5

    fig, axs = plt.subplots(1, 2, figsize=(definitions.FIGWIDTH, definitions.FIGWIDTH * aspect_performance), sharex=True, facecolor="white")

    cpu_hatch = "///"
    gpu_hatch = "\\"
    no_hatch = ""
    hatch_color = "#000000"
    color_grey = "#5B5B5B"
    color_fvm = "#CC6677"
    color_lbm = "#44AA99"
    color_cpu = "#0073FF"
    color_gpu = "#00FF95"

    axs[0].grid(True, linestyle="--")
    axs[0].bar(x[0], time_fvm, width=width, facecolor=color_fvm, edgecolor=color_fvm)
    axs[0].bar(x[1], time_lbm, width=width, facecolor=color_lbm, edgecolor=color_lbm)
    axs[0].set_yscale("log")
    axs[0].set_ylabel("Time-to-solution in s/rev")
    axs[0].set_ylim(1.0, 10000.0)
    axs[0].set_axisbelow(True)

    axs[1].grid(True, linestyle="--")

    # FVM: base colored bar + black hatch overlay
    axs[1].bar(x[0], energy_fvm, width=width, facecolor=color_fvm, edgecolor=color_fvm, lw=2.0, zorder=1)
    axs[1].bar(x[0], energy_fvm, width=width, facecolor="none", edgecolor=hatch_color, hatch=cpu_hatch, lw=0.0, zorder=2)

    # LBM GPU: base colored bar + black hatch overlay
    axs[1].bar(x[1], energy_lbm_gpu, width=width, facecolor=color_lbm, edgecolor=color_lbm, lw=2.0, zorder=1)
    axs[1].bar(x[1], energy_lbm_gpu, width=width, facecolor="none", edgecolor=hatch_color, hatch=gpu_hatch, lw=0.0, zorder=2)

    # LBM CPU stacked part: base colored bar + black hatch overlay
    if energy_lbm_cpu > 0.0:
        axs[1].bar(x[1], energy_lbm_cpu, width=width, bottom=energy_lbm_gpu, facecolor=color_lbm, edgecolor=color_lbm, lw=2.0, zorder=1)
        axs[1].bar(x[1], energy_lbm_cpu, width=width, bottom=energy_lbm_gpu, facecolor="none", edgecolor=hatch_color, hatch=cpu_hatch, lw=0.0, zorder=2)

    axs[1].set_yscale("log")
    axs[1].set_ylabel("Energy-to-solution in Wh/rev")
    axs[1].set_ylim(0.1, 1000.0)
    axs[1].set_axisbelow(True)

    for ax in axs:
        ax.set_xticks(x)
        ax.set_xticklabels(labels)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="white", hatch=cpu_hatch, edgecolor="black", label="CPU"),
        plt.Rectangle((0, 0), 1, 1, facecolor="white", hatch=gpu_hatch, edgecolor="black", label="GPU"),
    ]
    axs[1].legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout()
    fig.savefig(path_out, dpi=300, facecolor="white")
    plt.close(fig)
    return path_out