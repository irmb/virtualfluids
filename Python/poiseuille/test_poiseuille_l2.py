import unittest

import pyvista as pv
from norms import l2_norm
from poiseuille.analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille.simulation import run_simulation, grid_params, physical_params, sim_params
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices


class TestPoiseuilleFlow(unittest.TestCase):

    def test_poiseuille_flow(self):
        run_simulation()
        file_name = _get_output_file_name(sim_params)
        mesh = pv.read(file_name)
        indices = vertical_column_from_mesh(mesh)
        numerical_results = get_values_from_indices(mesh.get_array("Vx"), indices)
        heights = _get_heights_from_indices(mesh, indices)

        settings = _get_analytical_poiseuille_settings(grid_params, physical_params)
        analytical_results = poiseuille_at_heights(settings, heights)

        norm = l2_norm(analytical_results, numerical_results)
        print(f"L2 norm value: {norm}")
        self.assertLessEqual(norm, 1e-4)


def _get_analytical_poiseuille_settings(grid_params, physical_params):
    settings = PoiseuilleSettings()
    settings.length = grid_params.number_of_nodes_per_direction[0]
    settings.height = grid_params.number_of_nodes_per_direction[2]
    settings.viscosity = physical_params.lattice_viscosity
    settings.density = 1
    settings.force = 1e-6
    return settings


def _get_output_file_name(sim_params):
    file_name = f"output/mq/mq{sim_params.number_of_timesteps}/mq0_{sim_params.number_of_timesteps}.bin.vtu"
    return file_name


def _get_heights_from_indices(mesh, indices):
    return [mesh.points[index][2] for index in indices]

