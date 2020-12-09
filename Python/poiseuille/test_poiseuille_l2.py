import unittest

import pyvista as pv
from norms import l2_norm
from poiseuille.analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille.simulation import run_simulation, grid_params, physical_params, runtime_params
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices


class TestPoiseuilleFlow(unittest.TestCase):

    def test_poiseuille_flow(self):
        """
        WHEN comparing the simulation results to the analytical solution THEN the L2-Norm should be less than 1e-4
        """

        run_simulation()
        file_name_of_last_timestep = get_output_file_name(runtime_params)
        mesh_of_last_timestep = pv.read(file_name_of_last_timestep)
        column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
        numerical_results_from_single_column = get_values_from_indices(mesh_of_last_timestep.get_array("Vx"), column_indices)

        heights = get_heights_from_indices(mesh_of_last_timestep, column_indices)
        settings = get_analytical_poiseuille_settings(grid_params, physical_params)
        analytical_results = poiseuille_at_heights(settings, heights)

        l2_norm_result = l2_norm(analytical_results, numerical_results_from_single_column)

        max_acceptable_error = 1e-4
        self.assertLessEqual(l2_norm_result, max_acceptable_error)
        print(f"The L2-Norm is: {l2_norm_result}")


def get_analytical_poiseuille_settings(grid_params, physical_params):
    settings = PoiseuilleSettings()
    settings.length = grid_params.number_of_nodes_per_direction[0]
    settings.height = grid_params.number_of_nodes_per_direction[2]
    settings.viscosity = physical_params.lattice_viscosity
    settings.density = 1
    settings.force = 1e-6
    return settings


def get_output_file_name(rumtime_params):
    timesteps = rumtime_params.number_of_timesteps
    file_name = f"output/mq/mq{timesteps}/mq0_{timesteps}.bin.vtu"
    return file_name


def get_heights_from_indices(mesh, indices):
    return [mesh.points[index][2] for index in indices]

