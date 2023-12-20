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

! \author Sven Marcus, Henry Korb
=======================================================================================
"""
import unittest
from typing import List

from pyfluids import cpu
from acousticscaling import OneDirectionalAcousticScaling


class OneDirectionalAcousticScalingTest(unittest.TestCase):

    def setUp(self) -> None:
        self.grid_params = self.make_grid_params()
        self.physical_params = self.make_physical_params()
        self.runtime_params = self.make_runtime_params()
        self.kernel = self.make_kernel()

        self.sut = OneDirectionalAcousticScaling(self.grid_params, self.physical_params, self.runtime_params,
                                                 self.kernel)

    def test_given_sim_parameters__when_scaling_level_zero__should_return_equal_sim_parameters(self):
        factor = 1
        actual_params = self.sut.configuration_for_scale_level(0)
        actual_grid_params = actual_params[0]
        actual_physical_params = actual_params[factor]
        actual_runtime_params = actual_params[2]
        actual_kernel = actual_params[3]

        self.assert_parameters_scaled_by_factor(actual_grid_params, actual_kernel,
                                                actual_physical_params, actual_runtime_params, factor)

    def test_given_sim_parameters__when_scaling_level_one__should_return_sim_parameters_scaled_by_two(self):
        actual_params = self.sut.configuration_for_scale_level(1)
        actual_grid_params = actual_params[0]
        actual_physical_params = actual_params[1]
        actual_runtime_params = actual_params[2]
        actual_kernel = actual_params[3]

        self.assert_parameters_scaled_by_factor(actual_grid_params, actual_kernel,
                                                actual_physical_params, actual_runtime_params, 2)

    def assert_parameters_scaled_by_factor(self, actual_grid_params, actual_kernel,
                                           actual_physical_params, actual_runtime_params, factor):
        self.assert_grid_params_scaled_by_factor(actual_grid_params, factor=factor)
        self.assert_physical_params_scaled_by_factor(actual_physical_params, factor=factor)
        self.assert_runtime_params_scaled_by_factor(actual_runtime_params, factor=factor)
        self.assert_kernel_forcing_scaled_by_factor(actual_kernel, factor=factor)

    def assert_grid_params_scaled_by_factor(self, actual_grid_params: GridParameters, factor: int):
        expected_nodes_per_direction = self.scaled_list(self.grid_params.number_of_nodes_per_direction, factor)
        expected_blocks_per_direction = self.scaled_list(self.grid_params.blocks_per_direction, factor)
        expected_node_distance = self.grid_params.node_distance / factor
        self.assertEqual(expected_node_distance, actual_grid_params.node_distance)
        self.assertEqual(expected_nodes_per_direction, actual_grid_params.number_of_nodes_per_direction)
        self.assertEqual(expected_blocks_per_direction, actual_grid_params.blocks_per_direction)
        self.assertEqual(self.grid_params.reference_direction_index, actual_grid_params.reference_direction_index)
        self.assertEqual(self.grid_params.periodic_boundary_in_x1, actual_grid_params.periodic_boundary_in_x1)
        self.assertEqual(self.grid_params.periodic_boundary_in_x2, actual_grid_params.periodic_boundary_in_x2)
        self.assertEqual(self.grid_params.periodic_boundary_in_x3, actual_grid_params.periodic_boundary_in_x3)

    def assert_physical_params_scaled_by_factor(self, actual_params: cpu.parameters.PhysicalParameters, factor: int):
        self.assertEqual(self.physical_params.lattice_viscosity * factor, actual_params.lattice_viscosity)
        self.assertEqual(self.physical_params.bulk_viscosity_factor, actual_params.bulk_viscosity_factor)

    def assert_runtime_params_scaled_by_factor(self, actual_params: cpu.parameters.RuntimeParameters, factor: int):
        self.assertEqual(self.runtime_params.number_of_timesteps * factor, actual_params.number_of_timesteps)
        self.assertEqual(self.runtime_params.number_of_threads, actual_params.number_of_threads)
        self.assertEqual(self.runtime_params.timestep_log_interval, actual_params.timestep_log_interval)

    def assert_kernel_forcing_scaled_by_factor(self, actual_kernel: cpu.kernel.LBMKernel, factor: int):
        self.assertEqual(self.kernel.type, actual_kernel.type)
        self.assertEqual(self.kernel.use_forcing, actual_kernel.cpu.parameters.use_forcing)
        self.assertAlmostEqual(self.kernel.forcing_in_x1 / factor, actual_kernel.forcing_in_x1)
        self.assertAlmostEqual(self.kernel.forcing_in_x2, actual_kernel.forcing_in_x2)
        self.assertAlmostEqual(self.kernel.forcing_in_x3, actual_kernel.forcing_in_x3)

    @staticmethod
    def scaled_list(list_to_scale: List[int], factor: int) -> List[int]:
        return [list_to_scale[0], list_to_scale[1], list_to_scale[2] * factor]

    @staticmethod
    def make_kernel():
        kernel = cpu.kernel.LBMKernel(cpu.kernel.KernelType.CompressibleCumulantFourthOrderViscosity)
        kernel.use_forcing = True
        kernel.forcing_in_x1 = 5e-10
        return kernel

    @staticmethod
    def make_runtime_params():
        runtime_params = cpu.parameters.RuntimeParameters()
        runtime_params.number_of_threads = 4
        runtime_params.number_of_timesteps = 4_000_000
        runtime_params.timestep_log_interval = 1_000_000
        return runtime_params

    @staticmethod
    def make_physical_params():
        physical_params = cpu.parameters.PhysicalParameters()
        physical_params.lattice_viscosity = 1e-4
        return physical_params

    @staticmethod
    def make_grid_params():
        grid_params = cpu.parameters.GridParameters()
        grid_params.node_distance = 1
        grid_params.number_of_nodes_per_direction = [1, 1, 16]
        grid_params.blocks_per_direction = [1, 1, 16]
        grid_params.periodic_boundary_in_x1 = True
        grid_params.periodic_boundary_in_x2 = True

        return grid_params


if __name__ == '__main__':
    unittest.main()
