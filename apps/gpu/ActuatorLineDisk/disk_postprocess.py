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

import argparse
import ast
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

DEFAULT_CONFIG_NAME = "apps/gpu/ActuatorLineDisk/actuatorlinedisk.cfg"
DEFAULT_FIGWIDTH = 6.0
PLANE_FILE_PATTERN = re.compile(r"^zPlane_bin_lev_(?P<level>\d+)_ID_\d+_Part_\d+_t_(?P<timestep>\d+)\.vtk\.bin\.vtu$")

try:
    import vf_postproc.definitions as vf_definitions

    FIGWIDTH = float(vf_definitions.FIGWIDTH)
except ModuleNotFoundError:
    FIGWIDTH = DEFAULT_FIGWIDTH

try:
    plt.style.use(["vf_postproc.vf_style"])
except OSError:
    pass


def _load_config(path: Path) -> dict[str, str]:
    config: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key and value:
            config[key] = value
    return config


def _get_nested_float_array(config: dict[str, str], key: str) -> np.ndarray:
    return np.asarray(ast.literal_eval(f"[{config[key]}]"), dtype=float)


def _get_float_array(config: dict[str, str], key: str) -> np.ndarray:
    return np.atleast_1d(np.asarray(ast.literal_eval(config[key]), dtype=float))


def _get_thrust_coefficients(config: dict[str, str]) -> np.ndarray:
    if "THRUST_COEFFICIENTS" in config:
        return _get_float_array(config, "THRUST_COEFFICIENTS")
    if "THRUST_COEFFICIENT" in config:
        return _get_float_array(config, "THRUST_COEFFICIENT")
    raise KeyError("Missing THRUST_COEFFICIENT or THRUST_COEFFICIENTS in config.")


def _resolve_input_path(path: Path, config_path: Path) -> Path:
    if path.is_absolute():
        return path

    candidates = [
        Path.cwd() / path,
        config_path.parent / path,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    return candidates[0]


def _format_case_tag(thrust_coefficient: float) -> str:
    return f"ct_{thrust_coefficient:.2f}".replace(".", "p")


def _format_ct_label(thrust_coefficient: float) -> str:
    return f"{thrust_coefficient:.2f}".rstrip("0").rstrip(".")


def _get_disk_center(disk_center: np.ndarray) -> np.ndarray:
    center = np.asarray(disk_center, dtype=float).reshape(-1)
    if center.size < 3:
        raise ValueError("POSITION_DISK_CENTER must contain x, y and z coordinates.")
    return center[:3]


def _get_case_path(path_planes: Path, thrust_coefficient: float) -> Path:
    case_tag = _format_case_tag(thrust_coefficient)
    direct_case_path = path_planes / case_tag
    if (direct_case_path / "zplane").exists():
        return direct_case_path
    if path_planes.name == case_tag and (path_planes / "zplane").exists():
        return path_planes
    return direct_case_path


def _get_latest_plane_file(path_planes_case: Path) -> Path:
    plane_files = []
    for candidate in path_planes_case.glob("zPlane_bin_lev_*_ID_*_Part_*_t_*.vtk.bin.vtu"):
        match = PLANE_FILE_PATTERN.match(candidate.name)
        if match is not None:
            plane_files.append((int(match.group("level")), int(match.group("timestep")), candidate))

    if not plane_files:
        raise FileNotFoundError(f"No z-plane probe files found in {path_planes_case}")

    _, _, plane_file = max(plane_files)
    return plane_file


def load_velocities(path_planes, thrust_coefficients, velocity_inlet, disk_center, diameter):
    path_planes = Path(path_planes)
    disk_center = _get_disk_center(disk_center)
    thrust_coefficients = np.asarray(thrust_coefficients, dtype=float)

    dict_vx = {"lbm": {}, "analytic_disk": {}, "analytic_far": {}}
    induction = 0.5 * (1.0 - np.sqrt(1.0 - thrust_coefficients))
    velocities_disk = 1.0 - induction
    velocities_far = 1.0 - 2.0 * induction

    for index, ct in enumerate(thrust_coefficients):
        path_planes_case = _get_case_path(path_planes, float(ct)) / "zplane"
        plane_file = _get_latest_plane_file(path_planes_case)
        plane = pv.read(plane_file)

        points = np.asarray(plane.points, dtype=float)
        x = points[:, 0]
        y = points[:, 1]
        vx_mean = np.asarray(plane.point_data["vx_mean"], dtype=float)

        coord_y = np.round(y, 10)
        y_unique = np.unique(coord_y)
        y_line = y_unique[np.argmin(np.abs(y_unique - disk_center[1]))]
        mask_y = np.isclose(coord_y, y_line)

        x_raw = x[mask_y] / float(diameter)
        vx_raw = vx_mean[mask_y] / float(velocity_inlet)

        x_unique = np.unique(x_raw)
        vx_profile = np.array([vx_raw[x_raw == x_value].mean() for x_value in x_unique], dtype=float)

        idx = np.argsort(x_unique)
        x_profile = x_unique[idx]
        vx_profile = vx_profile[idx]

        key = _format_ct_label(float(ct))
        dict_vx["lbm"][key] = np.column_stack((x_profile, vx_profile))
        dict_vx["analytic_disk"][key] = float(velocities_disk[index])
        dict_vx["analytic_far"][key] = float(velocities_far[index])

    return dict_vx


def plot_velocities(dict_vx, path_out, disk_center, diameter=1.0):
    disk_x = _get_disk_center(disk_center)[0] / float(diameter)
    path_out = Path(path_out)
    path_out.parent.mkdir(parents=True, exist_ok=True)

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

    xlim = [-10, 10]
    ylim = [0.65, 1.05]
    xticks = [-10, 0, 10]

    aspect_velocities = 0.2
    aspect_subplot = 0.4

    xlegend, ylegend = 0.80, 0.50

    margin_left = 0.1
    margin_right = 0.85
    margin_bottom = 0.22
    margin_top = 0.82
    wspace = -0.2
    hspace = 0.0

    ct_keys = list(dict_vx["lbm"].keys())
    if not ct_keys:
        raise ValueError("No velocity data available for plotting.")

    fig, axes = plt.subplots(
        nrows=1,
        ncols=len(ct_keys),
        sharex=True,
        sharey=True,
        figsize=(FIGWIDTH, FIGWIDTH * aspect_velocities),
        facecolor="white",
    )
    axes = np.atleast_1d(axes)

    fig.subplots_adjust(
        left=margin_left,
        right=margin_right,
        bottom=margin_bottom,
        top=margin_top,
        wspace=wspace,
        hspace=hspace,
    )

    for index, (ax, ct_key) in enumerate(zip(axes, ct_keys)):
        values = np.asarray(dict_vx["lbm"][ct_key], dtype=float)
        x = values[:, 0]
        u = values[:, 1]

        u_disk = float(dict_vx["analytic_disk"][ct_key])
        u_far = float(dict_vx["analytic_far"][ct_key])

        ax.axhline(
            u_disk,
            color=color_disk_velocity,
            linestyle=linestyle_disk_velocity,
            linewidth=linewidth_disk_velocity,
            label=label_disk if index == 0 else None,
        )
        ax.axhline(
            u_far,
            color=color_disk_far,
            linestyle=linestyle_disk_far,
            linewidth=linewidth_disk_far,
            label=label_far if index == 0 else None,
        )
        ax.plot(
            x - disk_x,
            u,
            color=color_lbm,
            linestyle=linestyle_lbm,
            linewidth=linewidth_lbm,
            label=label_lbm if index == 0 else None,
        )
        ax.axvline(
            0.0,
            color=color_disk_position,
            linestyle=linestyle_disk_position,
            linewidth=linewidth_disk_position,
        )

        ax.set_title(rf"$C_T = {ct_key}$")
        ax.text(-0.5 * disk_x, 1.045, "disk position", ha="center", va="top", color=color_disk_position)
        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_ylim(ylim)
        ax.set_xlabel(r"$x/D$")
        ax.set_box_aspect(aspect_subplot)

        if index > 0:
            ax.tick_params(axis="y", left=False, labelleft=False)

    axes[0].set_ylabel(r"$u/u_\infty$")

    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="center left", ncol=1, frameon=False, bbox_to_anchor=(xlegend, ylegend))

    fig.savefig(path_out, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Post-process ActuatorLineDisk probe-plane output.")
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG_NAME,
        help="Path to the ActuatorLineDisk configuration file. Defaults to actuatorlinedisk.cfg.",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    config_path = _resolve_input_path(Path(args.config), script_dir)
    config = _load_config(config_path)

    thrust_coefficients = _get_thrust_coefficients(config)
    velocity_inlet = float(config["VELOCITY_INLET"])
    disk_diameter = float(config["DISK_DIAMETER"])
    disk_center = _get_nested_float_array(config, "POSITION_DISK_CENTER") * disk_diameter

    output_case_path = _resolve_input_path(Path(config["Path"]), config_path)
    output_root = output_case_path.parent
    plot_path = output_root / "post" / "velocities.png"

    dict_vx = load_velocities(output_root, thrust_coefficients, velocity_inlet, disk_center, disk_diameter)
    plot_velocities(dict_vx, plot_path, disk_center, disk_diameter)

    print(f"Wrote validation plot to {plot_path}")


if __name__ == "__main__":
    main()
