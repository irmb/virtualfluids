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

=======================================================================================
"""
import os

from pyfluids import cpu
from pymuparser import Parser


def get_max_length(number_of_nodes_per_direction, delta_x):
    return (number_of_nodes_per_direction[0] * delta_x,
            number_of_nodes_per_direction[1] * delta_x,
            number_of_nodes_per_direction[2] * delta_x)


physical_params = cpu.parameters.PhysicalParameters()
physical_params.lattice_viscosity = 0.005

grid_params = cpu.parameters.GridParameters()
grid_params.number_of_nodes_per_direction = [200, 120, 120]
grid_params.blocks_per_direction = [2, 2, 2]
grid_params.node_distance = 0.125
grid_params.periodic_boundary_in_x1 = False
grid_params.periodic_boundary_in_x2 = True
grid_params.periodic_boundary_in_x3 = True

runtime_params = cpu.parameters.RuntimeParameters()
runtime_params.timestep_log_interval = 1000
runtime_params.number_of_timesteps = 50000
runtime_params.number_of_threads = int(os.environ.get("OMP_NUM_THREADS", 4))


def run_simulation(physical_parameters=physical_params, grid_parameters=grid_params,
                   runtime_parameters=runtime_params):
    wall_thickness = 3 * grid_parameters.node_distance

    min_x, min_y, min_z = 0, 0, 0
    max_x, max_y, max_z = get_max_length(grid_parameters.number_of_nodes_per_direction, grid_parameters.node_distance)

    bottom_wall = cpu.geometry.GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, min_z, max_x + wall_thickness,
                             max_y + wall_thickness, min_z - wall_thickness)

    top_wall = cpu.geometry.GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, max_z, max_x + wall_thickness,
                          max_y + wall_thickness,
                          max_z + wall_thickness)

    left_wall = cpu.geometry.GbCuboid3D(min_x - wall_thickness, min_y, min_z - wall_thickness, max_x + wall_thickness,
                           min_y - wall_thickness,
                           max_z + wall_thickness)

    right_wall = cpu.geometry.GbCuboid3D(min_x - wall_thickness, max_y, min_z - wall_thickness, max_x + wall_thickness,
                            max_y + wall_thickness, max_z + wall_thickness)

    obstacle = cpu.geometry.GbCuboid3D(7, 7, 7, 8, 8, 8)

    velocity_boundary = cpu.geometry.GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, min_z - wall_thickness, min_x,
                                   max_y + wall_thickness, max_z + wall_thickness)

    outflow_boundary = cpu.geometry.GbCuboid3D(max_x, min_y - wall_thickness, min_z - wall_thickness, max_x + wall_thickness,
                                  max_y + wall_thickness, max_z + wall_thickness)

    no_slip_bc = cpu.boundaryconditions.NoSlipBoundaryCondition()

    outflow_bc = cpu.boundaryconditions.DensityBoundaryCondition()

    velocity_function = Parser()
    velocity_function.define_constant("u", 0.07)
    velocity_function.expression = "u"
    velocity_bc = cpu.boundaryconditions.VelocityBoundaryCondition(True, False, False, velocity_function, 0, -10)

    kernel = cpu.kernel.LBMKernel(cpu.kernel.KernelType.CompressibleCumulantFourthOrderViscosity)
    # kernel.use_forcing = True
    # kernel.forcing_in_x1 = 3e-6

    writer = cpu.writer.Writer()
    writer.output_path = "./output"
    writer.output_format = cpu.writer.OutputFormat.BINARY

    simulation = cpu.Simulation()
    simulation.set_writer(writer)

    simulation.set_physical_parameters(physical_parameters)
    simulation.set_grid_parameters(grid_parameters)
    simulation.set_runtime_parameters(runtime_parameters)
    simulation.set_kernel_config(kernel)

    # simulation.add_object(bottom_wall, no_slip_bc, 1, "/geo/bottomWall")
    # simulation.add_object(top_wall, no_slip_bc, 1, "/geo/topWall")
    # simulation.add_object(left_wall, no_slip_bc, 1, "/geo/leftWall")
    # simulation.add_object(right_wall, no_slip_bc, 1, "/geo/rightWall")

    simulation.add_object(obstacle, no_slip_bc, 1, "/geo/obstacle")

    simulation.add_object(outflow_boundary, outflow_bc, 1, "/geo/outflow")
    simulation.add_object(velocity_boundary, velocity_bc, 1, "/geo/velocityBoundary")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
