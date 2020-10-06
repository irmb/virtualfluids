import unittest

import pyvista as pv
from norms import l2_norm
from poiseuille_analytical import poiseuille_at_heights, PoiseuilleSettings
from poiseuille_flow import simulate_poiseuille_flow
from pyfluids.parameters import PhysicalParameters, GridParameters, SimulationParameters
from vtk_utilities import vertical_column_from_mesh, get_values_from_indices


class TestPoiseuilleFlow(unittest.TestCase):

    def test_poiseuille_flow(self):
        physical_params = PhysicalParameters()
        physical_params.lattice_viscosity = 0.005

        grid_params = GridParameters()
        grid_params.delta_x = 1
        grid_params.number_of_nodes_per_direction = [2, 2, 10]
        grid_params.blocks_per_direction = [1, 1, 1]
        grid_params.periodic_boundary_in_x1 = True
        grid_params.periodic_boundary_in_x2 = True

        sim_params = SimulationParameters()
        sim_params.number_of_threads = 4
        sim_params.number_of_timesteps = 10000
        sim_params.timestep_log_interval = 1000

        simulate_poiseuille_flow(physical_params, grid_params, sim_params)
        file_name = f"output/mq/mq{sim_params.number_of_timesteps}/mq0_{sim_params.number_of_timesteps}.bin.vtu"
        mesh = pv.read(file_name)
        indices = vertical_column_from_mesh(mesh)
        numerical_results = get_values_from_indices(mesh.get_array("Vx"), indices)
        heights = [mesh.points[index][2] for index in indices]

        settings = PoiseuilleSettings()
        settings.length = grid_params.number_of_nodes_per_direction[0]
        settings.height = grid_params.number_of_nodes_per_direction[2]
        settings.viscosity = physical_params.lattice_viscosity
        settings.density = 1
        settings.force = 1e-6

        analytical_results = poiseuille_at_heights(settings, heights)

        norm = l2_norm(analytical_results, numerical_results)
        print(f"L2 norm value: {norm}")
        self.assertLessEqual(norm, 1e-4)
