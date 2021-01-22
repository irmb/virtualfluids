import shutil
import unittest

import pyvista as pv
from pyfluids.parameters import GridParameters, PhysicalParameters, RuntimeParameters

from norms import root_mean_squared_error
from poiseuille.analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille.simulation import run_simulation
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices


class TestPoiseuilleFlow(unittest.TestCase):

    def test_poiseuille_flow(self):
        """
        WHEN comparing the simulation results to the analytical solution THEN the L2-Norm should be less than 1e-4
        """
        self.skipTest("Skipping test! This test is not implemented correctly")

        physical_params = PhysicalParameters()
        physical_params.lattice_viscosity = 0.0005

        runtime_params = RuntimeParameters()
        runtime_params.number_of_threads = 4
        runtime_params.timestep_log_interval = 1000

        runtime_params.number_of_timesteps = 10000
        channel_height = 10
        nodes_in_column = 8
        grid_params = create_grid_params_with_nodes_in_column(nodes_in_column,
                                                              delta_x=channel_height / nodes_in_column)
        l2_norm_result_100 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, channel_height)

        runtime_params.number_of_timesteps *= 2
        physical_params.lattice_viscosity *= 2
        nodes_in_column *= 2
        grid_params = create_grid_params_with_nodes_in_column(nodes_in_column,
                                                              delta_x=channel_height / nodes_in_column)
        l2_norm_result_200 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, channel_height)

        runtime_params.number_of_timesteps *= 2
        physical_params.lattice_viscosity *= 2
        nodes_in_column *= 2
        grid_params = create_grid_params_with_nodes_in_column(nodes_in_column,
                                                              delta_x=channel_height / nodes_in_column)
        l2_norm_result_400 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, channel_height)


        self.assertTrue(l2_norm_result_200 <= l2_norm_result_100)
        self.assertTrue(l2_norm_result_400 <= l2_norm_result_200)


def get_delta_x(number_of_nodes, height):
    return height / number_of_nodes


def run_simulation_with_settings(grid_params, physical_params, runtime_params, output_folder):
    remove_existing_output_directory(output_folder)
    run_simulation(physical_params, grid_params, runtime_params)


def get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, channel_height):
    output_folder = "./output"
    run_simulation_with_settings(grid_params, physical_params, runtime_params, output_folder)
    heights = get_heights(output_folder, runtime_params)

    numerical_results = get_numerical_results(runtime_params, output_folder, heights)
    analytical_results = get_analytical_results(physical_params, heights, channel_height)

    return root_mean_squared_error(analytical_results, numerical_results)


def get_heights(output_folder, runtime_params):
    mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
    column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
    heights = get_heights_from_indices(mesh_of_last_timestep, column_indices)
    return heights


def get_numerical_results(runtime_params, output_folder, heights):
    mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
    velocities_in_x_direction = mesh_of_last_timestep.get_array("Vx")
    column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
    numerical_results = get_values_from_indices(velocities_in_x_direction, column_indices)

    return numerical_results


def calculate_analytical_results(physical_params, height_values, channel_height):
    settings = get_analytical_poiseuille_settings(channel_height, physical_params)
    analytical_results = poiseuille_at_heights(settings, height_values)
    return analytical_results


def get_analytical_results(physical_params, heights, channel_height):
    analytical_results = calculate_analytical_results(physical_params, heights, channel_height)
    return analytical_results


def get_mesh_for_last_timestep(output_folder, runtime_params):
    file_name_of_last_timestep = get_output_file_name(output_folder, runtime_params)
    mesh_of_last_timestep = pv.read(file_name_of_last_timestep)
    return mesh_of_last_timestep


def remove_existing_output_directory(output_dir):
    shutil.rmtree(output_dir, ignore_errors=True)


def get_analytical_poiseuille_settings(height, physical_params):
    settings = PoiseuilleSettings()
    settings.height = height
    settings.viscosity = physical_params.lattice_viscosity
    settings.density = 1
    settings.force = 1e-6
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
    grid_params.number_of_nodes_per_direction = [8, 8, nodes_in_column]
    grid_params.blocks_per_direction = [1, 1, 4]
    grid_params.periodic_boundary_in_x1 = True
    grid_params.periodic_boundary_in_x2 = True
    grid_params.periodic_boundary_in_x3 = False

    print(f"GridParameters.node_distance = {grid_params.node_distance}")
    print(f"GridParameters.number_of_nodes_per_direction = {grid_params.number_of_nodes_per_direction}")

    return grid_params
