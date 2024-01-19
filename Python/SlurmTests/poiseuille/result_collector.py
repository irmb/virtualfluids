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
# from typing import Collection, List

# import pyvista as pv
# from poiseuille.analytical import PoiseuilleSettings, poiseuille_at_heights
# from vtk_utilities import vertical_column_from_mesh, get_values_from_indices
# from SlurmTests.poiseuille.settings import Scaling


# def get_output_file_name(output_folder, runtime_params):
#     timesteps = runtime_params.number_of_timesteps
#     file_name = f"{output_folder}/mq/mq{timesteps}/mq0_{timesteps}.bin.vtu"

#     return file_name


# def get_mesh_for_last_timestep(output_folder, runtime_params):
#     file_name_of_last_timestep = get_output_file_name(output_folder, runtime_params)
#     mesh_of_last_timestep = pv.read(file_name_of_last_timestep)
#     return mesh_of_last_timestep


# def get_heights_from_indices(mesh, indices):
#     return [mesh.points[index][2] for index in indices]


# def get_heights(output_folder, runtime_params):
#     mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
#     column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
#     heights = get_heights_from_indices(mesh_of_last_timestep, column_indices)
#     return heights


# def get_numerical_results(runtime_params, output_folder):
#     mesh_of_last_timestep = get_mesh_for_last_timestep(output_folder, runtime_params)
#     velocities_in_x_direction = mesh_of_last_timestep.get_array("Vx")
#     column_indices = vertical_column_from_mesh(mesh_of_last_timestep)
#     numerical_results = get_values_from_indices(velocities_in_x_direction, column_indices)

#     return numerical_results


# def get_analytical_results(grid_params, physical_params, kernel, height_values):
#     channel_height = grid_params.number_of_nodes_per_direction[2]
#     settings = get_analytical_poiseuille_settings(channel_height, physical_params, kernel)
#     max_grid_height = channel_height * grid_params.node_distance
#     adjusted_height_values = [value / max_grid_height * channel_height for value in height_values]
#     analytical_results = poiseuille_at_heights(settings, adjusted_height_values)
#     return analytical_results


# def get_analytical_poiseuille_settings(height, physical_params, kernel):
#     settings = PoiseuilleSettings()
#     settings.height = height
#     settings.viscosity = physical_params.lattice_viscosity
#     settings.density = 1
#     settings.force = kernel.forcing_in_x1

#     return settings


# def collect_results() -> (List[List[float]], List[List[float]]):
#     analytical_results = []
#     numerical_results = []

#     for simulation_run in range(0, 3):
#         output_folder = f"output-{simulation_run}"
#         grid_params, physical_params, runtime_params, kernel = Scaling.configuration_for_scale_level(simulation_run)
#         heights = get_heights(output_folder, runtime_params)
#         analytical_results.append(
#             get_analytical_results(grid_params, physical_params, kernel, heights))
#         numerical_results.append(get_numerical_results(runtime_params, output_folder))

#     return analytical_results, numerical_results
