# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
# %%
"""
This script runs the Heated Cube app for different Rayleigh numbers and compares to reference data from 
https://www.sciencedirect.com/science/article/pii/S0017931000000375
"""
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
from pathlib import Path
import subprocess
from dataclasses import dataclass

# %%
thermal_expansion_velocity = 1
t_hot = 0.5
t_cold = -0.5
CMAP = "coolwarm"
FIGSIZE = 503 / 72.27 * 0.45

rayleigh_numbers = [1000, 10_000, 100_000, 1_000_000]
diffusivity_lb = 0.001
nx = 64
t_outs = [0.1, 0.2, 0.2, 0.2]

#%%
THIS_DIR = Path(__file__).parent
MAIN_DIR = THIS_DIR.parent.parent
OUTPUT_DIR = MAIN_DIR/"apps/gpu/HeatedCube/output"
executable = MAIN_DIR/"build/bin//HeatedCube"

#%%
for rayleigh_number, t_out in zip(rayleigh_numbers, t_outs):
    file_name = OUTPUT_DIR/f"heated_cubeRa1e{int(np.log10(rayleigh_number)):d}.cfg"
    with open(file_name, "w") as f:
        text = [
                f"Path = ./output/Ra1e{int(np.log10(rayleigh_number)):d}/",
                "tStartOut = 0",
                f"tOut = {t_out:1.3f}",
                f"tEnd = {t_out * 3.01:1.3f}",
                f"nNodes = {nx}",
                f"Ra = {rayleigh_number}",
                f"diffusivityLB = {diffusivity_lb}",
            ]
        f.writelines("\n".join(text))
    subprocess.run([executable, file_name])
#%%

plt.rcParams.update(
    {
        # set x axis
        "xtick.direction": "in",
        "xtick.major.size": 3,
        "xtick.major.width": 0.5,
        "xtick.minor.size": 1.5,
        "xtick.minor.width": 0.5,
        "xtick.minor.visible": True,
        "xtick.top": True,
        # Set y axis
        "ytick.direction": "in",
        "ytick.major.size": 3,
        "ytick.major.width": 0.5,
        "ytick.minor.size": 1.5,
        "ytick.minor.width": 0.5,
        "ytick.minor.visible": True,
        "ytick.right": True,
        # Set line widths
        "axes.linewidth": 0.5,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.0,
        # Remove legend frame
        "legend.frameon": False,
        # Always save as 'tight'
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        # Spacing around subfigures
        "figure.constrained_layout.w_pad": 0.01,
        "figure.constrained_layout.h_pad": 0.02,
        "figure.constrained_layout.wspace": 0,
        "figure.constrained_layout.hspace": 0.05,
        "figure.dpi": 300,
        "figure.constrained_layout.use": True,
        # Use serif fonts
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "text.parse_math": True,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 10,
        "font.size": 10,
    }
)

def get_midplane_filename(time_step: int | str):
    return f"midPlane_bin_lev_0_ID_0_Part_1_t_{time_step}.vtk.bin.vtu"


def get_sideplane_filename(time_step: int | str):
    return f"sidePlane_bin_lev_0_ID_0_Part_1_t_{time_step}.vtk.bin.vtu"

@dataclass
class Reference:
    name: str
    u_max: float
    x: float
    z: float
    nusselt_number_mid: float
    nusselt_number_3d: float

    def get_files(self, file_number: int) -> tuple[pv.UnstructuredGrid, pv.UnstructuredGrid]:
        output_dir =OUTPUT_DIR / self.name
        mid_plane_files = list(output_dir.glob(get_midplane_filename("*")))
        side_plane_files = list(output_dir.glob(get_sideplane_filename("*")))
        if not mid_plane_files:
            print(f"no plane files found in {output_dir}")
        sorted_mid_plane_files = sorted(mid_plane_files, key=lambda x: int(x.name.split("_")[-1].split(".")[0]))
        sorted_side_plane_files = sorted(side_plane_files, key=lambda x: int(x.name.split("_")[-1].split(".")[0]))
        return pv.read(sorted_mid_plane_files[file_number]), pv.read(sorted_side_plane_files[file_number])

    def plot_comparison(self, file_number: int, ax: plt.Axes):
        mid_plane, side_plane = self.get_files(file_number)
        nx = int(np.sqrt(mid_plane.points.shape[0]))
        velocity = mid_plane.get_array("vx").reshape((nx, nx)) / thermal_expansion_velocity
        temperature = mid_plane.get_array("phi").reshape((nx, nx))
        max_index = np.argmax(velocity)

        ax.scatter(mid_plane.points[max_index, 0], mid_plane.points[max_index, 2], color="k", marker="x")
        ax.scatter(self.x, self.z, facecolor="none", edgecolor="black", marker="o")
        ax.set_title(self.name)
        return ax.imshow(temperature, extent=(-0.5, 0.5, -0.5, 0.5), origin="lower", cmap=CMAP)

    def compare(self, file_number: int):
        mid_plane, side_plane = self.get_files(file_number)
        delta_x = mid_plane.points[1, 0] - mid_plane.points[0, 0]
        velocity = mid_plane.get_array("vx") / thermal_expansion_velocity
        nx = int(np.sqrt(velocity.shape[0]))

        max_index = np.argmax(velocity)
        max_value = velocity[max_index]
        max_location = mid_plane.points[max_index]
        value_error = self.u_max / max_value - 1
        x_error = np.abs(max_location[0] - self.x) / delta_x
        z_error = np.abs(max_location[2] - self.z) / delta_x
        temperature = side_plane.get_array("phi").reshape((nx, nx)).T
        nusselt_number_mid = np.mean((temperature[nx // 2, :] - t_cold) / (0.5 * delta_x))
        nusselt_number_3d = np.mean((temperature - t_cold) / (0.5 * delta_x))
        nusselt_error_mid = nusselt_number_mid / self.nusselt_number_mid - 1
        nusselt_error_3d = nusselt_number_3d / self.nusselt_number_3d - 1
        print(
            self.name,
            f"vx: {value_error * 100: 1.3f}%, x: {x_error: 1.3f} dx, z: {z_error: 1.3f} dx Nu_mid: {nusselt_number_mid: 1.5f} / {self.nusselt_number_mid:1.5f} {nusselt_error_mid * 100:1.5}% Nu3d {nusselt_number_3d:1.5f} / {self.nusselt_number_3d:1.5f} {nusselt_error_3d * 100:1.5f}%",
        )
        return temperature


references = [
    Reference("Ra1e3", 3.54356, 0.0166, 0.3169, 1.087, 1.07),
      Reference("Ra1e4", 16.71986, 0.0196, 0.3250, 2.2505, 2.0542),
      Reference("Ra1e5", 43.0610, -0.1865, 0.3848, 4.6127, 4.3370),
      Reference("Ra1e6", 123.4777, -0.3133, 0.4366, 8.8771, 8.6407),
    #   Reference("Ra1e7", 383.8358, -0.3777, 0.4662, 16.5477, 16.3427)
]
#%%
for i in range(2):
    for reference in references:
        reference.compare(-i)
# %%
fig, axes = plt.subplots(2, 2, figsize=(FIGSIZE, FIGSIZE * 1.1), sharex="col", sharey="row", layout="constrained")
for reference, ax in zip(references, axes.flatten()):
    cbar = reference.plot_comparison(-3, ax)
fig.colorbar(cbar, label=r"$\theta$", ax=axes[0, :], orientation="horizontal", location="top")
for ax in axes.T[0]:
    ax.set_ylabel("$z/L$")
for ax in axes[-1]:
    ax.set_xlabel("$y/L$")
fig.savefig(OUTPUT_DIR/"comparison.png")
fig.savefig(OUTPUT_DIR/"comparison.pdf")

# %%
