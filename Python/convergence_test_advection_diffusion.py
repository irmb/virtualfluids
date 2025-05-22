# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder#%%
# %%
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import subprocess
#%%
MAIN_DIR = Path(__file__).parent.parent
GAUSSIAN_HILL_DIR = MAIN_DIR / "apps/gpu/GaussianHillOfConcentration"
OUTPUT_DIR =  GAUSSIAN_HILL_DIR / "output"
BINARY = MAIN_DIR / "build/bin/GaussianHillOfConcentration"
PLOT = True
COLORS = ["#d95f02", "#1b9e77", "#7570b3"]
MARKERS = ["x", "+"]
FIG_WIDTH = 503 / 72.27 * 0.5
CMAP = "Blues"
RUN = True

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
        "figure.constrained_layout.h_pad": 0.01,
        "figure.constrained_layout.wspace": 0,
        "figure.constrained_layout.hspace": 0,
        "figure.dpi": 300,
        "figure.constrained_layout.use": True,
        # Use serif fonts
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "text.parse_math": True,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.labelsize": 10,
        "legend.fontsize": 10,
        "font.size": 10,
    }
)

c_0 = 1
sigma_0 = 1.0
delta_t = 1.0
ns_per_sigma = np.array([2, 4, 8, 16])
peclet_numbers = [1, 10000]
diffusivities_lb = np.array([1e-3, 1e-5])
transports = [0.3, 3]


def write_config_file(n_per_sigma, peclet_number, diff_lb, transport):
    fname = OUTPUT_DIR / f"gaussian_hill_Pe{int(peclet_number)}_N{n_per_sigma}.cfg"
    text = "\n".join(
        [
            f"Path = output/Pe{int(peclet_number)}N{n_per_sigma}/",
            f"PecletNumber = {peclet_number}",
            f"NodesPerSigma0 = {n_per_sigma}",
            f"DiffusivityLB = {diff_lb}",
            f"Transports = {transport}",
        ]
    )
    with open(fname, "w") as f:
        f.writelines(text)
    return fname


# %%
if RUN:
    files = []
    for n_per_sigma in ns_per_sigma:
        for pe, diff, transport in zip(peclet_numbers, diffusivities_lb, transports):
            cfg = write_config_file(n_per_sigma, pe, diff, transport)
            files.append(cfg)

    for file in files:
        subprocess.run(args=[str(BINARY), str(file)], cwd=str(GAUSSIAN_HILL_DIR))


# %%
def read_time_from_name(name: str):
    return int(name.split("_")[-1].split(".")[0])


def analytical_solution(
    c_0: float, sigma_0: float, velocity: float, diffusivity: float, time: float, coords: np.ndarray
):
    scale = sigma_0**2 + 2 * diffusivity * time
    exponent = np.sum((coords - velocity * time) ** 2, axis=-1) / (2 * scale)
    return np.power(sigma_0 / np.sqrt(scale), 3.0) * c_0 * np.exp(-exponent)


def reshape(field: np.ndarray, nx: int):
    if len(field.shape) == 2:
        n_fields = field.shape[-1]
        return field.reshape((nx + 2, nx + 2, nx + 2, n_fields))[1:-1, 1:-1, 1:-1, :]
    return field.reshape((nx + 2, nx + 2, nx + 2))[1:-1, 1:-1, 1:-1]


# %%
rms_errors: dict[int, np.ndarray] = dict()
max_errors: dict[int, np.ndarray] = dict()
fig, ax = plt.subplots(
    len(peclet_numbers),
    len(ns_per_sigma),
    sharey=True,
    sharex=True,
    figsize=(FIG_WIDTH, FIG_WIDTH * 0.7),
    squeeze=False,
)

for i, (peclet_number, diffusivity_lb) in enumerate(zip(peclet_numbers, diffusivities_lb)):
    rms_errors[peclet_number] = np.zeros(len(ns_per_sigma))
    max_errors[peclet_number] = np.zeros(len(ns_per_sigma))
    for j, n_per_sigma in enumerate(ns_per_sigma):
        delta_x = sigma_0 / n_per_sigma
        diffusivity = diffusivity_lb * (delta_x**2) / delta_t
        output_dir = OUTPUT_DIR / f"Pe{peclet_number}N{n_per_sigma}"
        final_timestep = sorted(output_dir.glob("*.vtu"))[-1]
        n_t = read_time_from_name(final_timestep.name)
        velocity = peclet_number * diffusivity / sigma_0
        time = n_t * delta_t
        domain = pv.read(final_timestep)
        nx = int(np.cbrt(domain.points.shape[0] - 1)) - 2
        coords = reshape(domain.points[1:], nx)
        analytical = analytical_solution(c_0, sigma_0, velocity, diffusivity, time, coords + velocity * time / 2)
        numerical = reshape(domain.get_array("conc")[1:], nx)
        diff = analytical - numerical
        x_max = int(velocity * time / 2 / delta_x + nx // 2)
        colormap = ax[i, j].imshow(
            np.abs(diff[x_max]), vmax=0.1, vmin=1e-6, extent=[-9, 9, -9, 9], norm="log", origin="lower", cmap=CMAP
        )
        rms_errors[peclet_number][j] = np.sqrt(np.mean(diff**2))
        max_errors[peclet_number][j] = np.max(np.abs(analytical - numerical))
        mass_num = np.sum(numerical) * nx**3
        mass_an = np.sum(analytical) * nx**3
        print(mass_num / mass_an - 1)
for axis in ax.T[0]:
    axis.set_ylabel(r"$y/\sigma_0$")
for axis in ax[-1]:
    axis.set_xlabel(r"$x/\sigma_0$")
fig.colorbar(colormap, ax=ax, orientation="horizontal", location="top").set_label(r"$\Vert\Delta \theta\Vert$")
if PLOT:
    fig.savefig(OUTPUT_DIR / "differences.png")
    fig.savefig(OUTPUT_DIR / "differences.pdf")

# %%
for peclet_number in peclet_numbers:
    rms_fit = np.polynomial.polynomial.polyfit(np.log(ns_per_sigma), np.log(rms_errors[peclet_number]), 1)
    max_fit = np.polynomial.polynomial.polyfit(np.log(ns_per_sigma), np.log(max_errors[peclet_number]), 1)
    print(f"RMS order of convergence: {-rms_fit[1]:2.2f}")
    print(f"Max order of convergence: {-max_fit[1]:2.2f}")


# %%


def plot_error(
    cell_widths: np.ndarray,
    peclet_numbers: list[int],
    errors: dict[int, np.ndarray],
    ylabel: str,
    name: str,
):
    plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH * 0.6))
    x = np.linspace(3, 13, 2)
    plt.loglog(x, np.power(x, -2) * 0.9, "k-", label=r"$O(\Delta x^2)$")
    plt.loglog(x - 1, np.power(x - 1, -3) * 0.5, "k--", label=r"$O(\Delta x^3)$")
    for peclet_number, color, marker in zip(peclet_numbers, COLORS, MARKERS):
        plt.loglog(
            cell_widths,
            errors[peclet_number],
            linestyle="",
            marker=marker,
            color=color,
            label=r"$\mathrm{Pe} = " + str(peclet_number) + "$",
        )
    plt.xlabel("$N$")
    plt.ylabel(ylabel)
    plt.legend(bbox_to_anchor=(1, 0.5), loc="center left", columnspacing=1, handletextpad=0.4, fontsize=10)
    plt.minorticks_off()
    plt.xticks(ns_per_sigma, [f"${n}$" for n in ns_per_sigma])
    plt.savefig(OUTPUT_DIR / (name + ".png"))
    plt.savefig(OUTPUT_DIR / (name + ".pdf"))


# %%
if PLOT:
    plot_error(ns_per_sigma, peclet_numbers, rms_errors, r"$\mathrm{rms}(\Delta \theta)$", "rms_error")
    plot_error(ns_per_sigma, peclet_numbers, max_errors, r"$\max\Vert \Delta \theta\Vert$", "max_error")

# %%
