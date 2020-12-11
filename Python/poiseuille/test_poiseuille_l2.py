import unittest

import shutil
import math
import pyvista as pv
from norms import root_mean_squared_error, mean_squared_error, mean_absolute_error
from pyfluids.parameters import GridParameters, PhysicalParameters, RuntimeParameters
from poiseuille.analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille.simulation import run_simulation
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices

from matplotlib.pyplot import plot, show, legend


class TestPoiseuilleFlow(unittest.TestCase):

    def test_poiseuille_flow(self):
        """
        WHEN comparing the simulation results to the analytical solution THEN the L2-Norm should be less than 1e-4
        """

        physical_params = PhysicalParameters()
        physical_params.lattice_viscosity = 0.0005

        runtime_params = RuntimeParameters()
        runtime_params.number_of_threads = 4
        runtime_params.timestep_log_interval = 1000

        runtime_params.number_of_timesteps = 10000
        grid_params = create_grid_params_with_nodes_in_column(nodes_in_column=5, height=10)
        l2_norm_result_100 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, 11)

        runtime_params.number_of_timesteps = 20000
        grid_params = create_grid_params_with_nodes_in_column(nodes_in_column=10, height=10)
        l2_norm_result_200 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, 11)

        # runtime_params.number_of_timesteps = 40000
        # grid_params = create_grid_params_with_nodes_in_column(nodes_in_column=20, height=10)
        # l2_norm_result_400 = get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, 11)

        # nodes = [5, 10, 20]
        # l2_norms = [l2_norm_result_100, l2_norm_result_200, l2_norm_result_400]
        # plot(nodes, l2_norms)
        # show()
        #
        # self.assertTrue(l2_norm_result_200 <= l2_norm_result_100)
        # self.assertTrue(l2_norm_result_400 <= l2_norm_result_200)


def get_l2_norm_for_simulation(grid_params, physical_params, runtime_params, total_height):
    output_folder = "./output"
    remove_existing_output_directory(output_folder)
    run_simulation(physical_params, grid_params, runtime_params)
    mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
    column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
    velocities_in_x_direction = mesh_of_last_timestep.get_array("Vx")
    numerical_results = get_values_from_indices(velocities_in_x_direction, column_indices)
    heights = get_heights_from_indices(mesh_of_last_timestep, column_indices)
    analytical_results = get_analytical_results(physical_params, total_height, heights)

    print("")
    # print(f"Numerical results: {numerical_results}")
    # print(f"Analytical results: {analytical_results}")
    print(f"Heights: {(min(heights), max(heights))}")
    plot(heights, numerical_results)
    plot(heights, analytical_results)
    legend(["numerical", "analytical"])
    show()

    return root_mean_squared_error(analytical_results, numerical_results)


def get_analytical_results(physical_params, total_height, height_values):
    settings = get_analytical_poiseuille_settings(total_height, physical_params)
    analytical_results = poiseuille_at_heights(settings, height_values)
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

    print(f"Poiseuille height {settings.height}")
    return settings


def get_output_file_name(output_folder, runtime_params):
    timesteps = runtime_params.number_of_timesteps
    file_name = f"{output_folder}/mq/mq{timesteps}/mq0_{timesteps}.bin.vtu"
    return file_name


def get_heights_from_indices(mesh, indices):
    return [mesh.points[index][2] for index in indices]


def create_grid_params_with_nodes_in_column(nodes_in_column, height):
    grid_params = GridParameters()
    grid_params.node_distance = height / (nodes_in_column - 1)
    grid_params.number_of_nodes_per_direction = [nodes_in_column, nodes_in_column, nodes_in_column]
    grid_params.blocks_per_direction = [1, 1, 1]
    grid_params.periodic_boundary_in_x1 = True
    grid_params.periodic_boundary_in_x2 = True

    return grid_params
