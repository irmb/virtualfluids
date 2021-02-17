import os
import shutil
import unittest

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import GridParameters, PhysicalParameters, RuntimeParameters
from scipy import stats

from errors import normalized_l2_error
from poiseuille.analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille.simulation import run_simulation
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices


class TestPoiseuilleFlow(unittest.TestCase):
    node_distances = [1, 0.5, 0.25]
    number_of_nodes = [16, 32, 64]
    number_of_timesteps = [2_500_000, 5_000_000, 10_000_000]
    forcings = [1e-9, 5e-10, 2.5e-10]
    viscosities = [1e-3, 2e-3, 4e-3]

    def zipped_settings(self):
        return zip(self.node_distances,
                   self.number_of_nodes,
                   self.number_of_timesteps,
                   self.forcings,
                   self.viscosities)

    def test_poiseuille_flow(self):
        self.skipTest("This test is not implemented correctly yet")
        plt.ion()

        physical_params = PhysicalParameters()

        runtime_params = RuntimeParameters()
        runtime_params.number_of_threads = os.cpu_count()
        runtime_params.timestep_log_interval = 10000

        kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
        kernel.use_forcing = True

        normalized_l2_errors = []
        for delta_x, nodes, timesteps, forcing, viscosity in self.zipped_settings():
            physical_params.lattice_viscosity = viscosity
            runtime_params.number_of_timesteps = timesteps
            kernel.forcing_in_x1 = forcing

            grid_params = create_grid_params_with_nodes_in_column(nodes, delta_x)
            l2_error = get_l2_error_for_simulation(grid_params, physical_params, runtime_params, kernel)
            normalized_l2_errors.append(l2_error)

        nodes_as_log = [np.log10(node) for node in self.number_of_nodes]
        l2_norms_as_log = [np.log10(l2) for l2 in normalized_l2_errors]
        res = stats.linregress(nodes_as_log, l2_norms_as_log)

        plt.xscale("log")
        plt.yscale("log")
        plt.plot(self.number_of_nodes, [np.power(10, res.intercept + res.slope * node) for node in nodes_as_log], 'r-')
        plt.plot(self.number_of_nodes, normalized_l2_errors, "x:")
        plt.show()

        print(normalized_l2_errors)
        self.assertAlmostEqual(res.slope, -2, places=2)


def get_l2_error_for_simulation(grid_params, physical_params, runtime_params, kernel):
    output_folder = "./output"
    run_simulation_with_settings(grid_params, physical_params, runtime_params, kernel, output_folder)
    heights = get_heights(output_folder, runtime_params)

    numerical_results = get_numerical_results(runtime_params, output_folder)
    analytical_results = get_analytical_results(grid_params, physical_params, kernel, heights)

    plt.plot(heights, numerical_results)
    plt.plot(heights, analytical_results)
    plt.legend(["numerical", "analytical"])
    plt.show()

    return normalized_l2_error(analytical_results, numerical_results)


def run_simulation_with_settings(grid_params, physical_params, runtime_params, kernel, output_folder):
    shutil.rmtree(output_folder, ignore_errors=True)
    run_simulation(physical_params, grid_params, runtime_params, kernel)


def get_heights(output_folder, runtime_params):
    mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
    column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
    heights = get_heights_from_indices(mesh_of_last_timestep, column_indices)
    return heights


def get_numerical_results(runtime_params, output_folder):
    mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
    velocities_in_x_direction = mesh_of_last_timestep.get_array("Vx")
    column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
    numerical_results = get_values_from_indices(velocities_in_x_direction, column_indices)

    return numerical_results


def get_analytical_results(grid_params, physical_params, kernel, height_values):
    channel_height = grid_params.number_of_nodes_per_direction[2]
    settings = get_analytical_poiseuille_settings(channel_height, physical_params, kernel)
    max_grid_height = channel_height * grid_params.node_distance
    adjusted_height_values = [value / max_grid_height * channel_height for value in height_values]
    analytical_results = poiseuille_at_heights(settings, adjusted_height_values)
    return analytical_results


def get_mesh_for_last_timestep(output_folder, runtime_params):
    file_name_of_last_timestep = get_output_file_name(output_folder, runtime_params)
    mesh_of_last_timestep = pv.read(file_name_of_last_timestep)
    return mesh_of_last_timestep


def get_analytical_poiseuille_settings(height, physical_params, kernel):
    settings = PoiseuilleSettings()
    settings.height = height
    settings.viscosity = physical_params.lattice_viscosity
    settings.density = 1
    settings.force = kernel.forcing_in_x1

    return settings


def get_output_file_name(output_folder, runtime_params):
    timesteps = runtime_params.number_of_timesteps
    file_name = f"{output_folder}/mq/mq{timesteps}/mq0_{timesteps}.bin.vtu"

    return file_name


def get_heights_from_indices(mesh, indices):
    return [mesh.points[index][2] for index in indices]


def create_grid_params_with_nodes_in_column(nodes_in_column, delta_x):
    grid_params = GridParameters()
    grid_params.node_distance = delta_x
    grid_params.number_of_nodes_per_direction = [1, 1, nodes_in_column]
    grid_params.blocks_per_direction = [1, 1, 8]
    grid_params.periodic_boundary_in_x1 = True
    grid_params.periodic_boundary_in_x2 = True
    grid_params.periodic_boundary_in_x3 = False

    return grid_params
