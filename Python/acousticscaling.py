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
from pyfluids import cpu


class OneDirectionalAcousticScaling:

    def __init__(self, grid_parameters: cpu.parameters.GridParameters,
                 physical_parameters: cpu.parameters.PhysicalParameters,
                 runtime_parameters: cpu.parameters.RuntimeParameters,
                 kernel: cpu.kernel.LBMKernel):
        self._grid_params = grid_parameters
        self._physical_params = physical_parameters
        self._runtime_params = runtime_parameters
        self._kernel = kernel

    def configuration_for_scale_level(self, level: int = 1) -> tuple[cpu.parameters.GridParameters,
                                                                cpu.parameters.PhysicalParameters,
                                                                cpu.parameters.RuntimeParameters,
                                                                cpu.kernel.LBMKernel]:
        if level < 0:
            raise ValueError("level must be >= 0")

        grid_params = self.clone_grid_params_for_level(level)
        physical_params = self.clone_physical_parameters(level)
        runtime_params = self.clone_runtime_params_for_level(level)
        kernel = self.clone_kernel_for_level(level)

        return grid_params, physical_params, runtime_params, kernel

    def clone_grid_params_for_level(self, level) -> cpu.parameters.GridParameters:
        grid_params = cpu.parameters.GridParameters()
        grid_params.reference_direction_index = self._grid_params.reference_direction_index
        grid_params.periodic_boundary_in_x1 = self._grid_params.periodic_boundary_in_x1
        grid_params.periodic_boundary_in_x2 = self._grid_params.periodic_boundary_in_x2
        grid_params.periodic_boundary_in_x3 = self._grid_params.periodic_boundary_in_x3

        grid_params.number_of_nodes_per_direction = list(self._grid_params.number_of_nodes_per_direction)
        grid_params.blocks_per_direction = list(self._grid_params.blocks_per_direction)
        grid_params.node_distance = self._grid_params.node_distance

        if level > 0:
            grid_params.node_distance /= (level * 2)
            grid_params.number_of_nodes_per_direction = [grid_params.number_of_nodes_per_direction[0],
                                                         grid_params.number_of_nodes_per_direction[1],
                                                         grid_params.number_of_nodes_per_direction[2] * (level * 2)]

            grid_params.blocks_per_direction = [grid_params.blocks_per_direction[0],
                                                grid_params.blocks_per_direction[1],
                                                grid_params.blocks_per_direction[2] * (level * 2)]

        return grid_params

    def clone_physical_parameters(self, level):
        physical_params = cpu.parameters.PhysicalParameters()
        physical_params.lattice_viscosity = self._physical_params.lattice_viscosity

        if level > 0:
            physical_params.lattice_viscosity *= (level * 2)

        return physical_params

    def clone_runtime_params_for_level(self, level):
        runtime_params = cpu.parameters.RuntimeParameters()
        runtime_params.number_of_timesteps = self._runtime_params.number_of_timesteps
        runtime_params.number_of_threads = self._runtime_params.number_of_threads
        runtime_params.timestep_log_interval = self._runtime_params.timestep_log_interval

        if level > 0:
            runtime_params.number_of_timesteps *= (level * 2)

        return runtime_params

    def clone_kernel_for_level(self, level):
        kernel = cpu.kernel.LBMKernel(self._kernel.type)
        kernel.use_forcing = self._kernel.use_forcing
        kernel.forcing_in_x1 = self._kernel.forcing_in_x1
        kernel.forcing_in_x2 = self._kernel.forcing_in_x2
        kernel.forcing_in_x3 = self._kernel.forcing_in_x3

        if level > 0:
            kernel.forcing_in_x1 /= (level * 2)
            kernel.forcing_in_x2 /= (level * 2)
            kernel.forcing_in_x3 /= (level * 2)

        return kernel
