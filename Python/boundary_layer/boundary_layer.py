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

# %%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger, parallel

# %%
sim_name = "ABL"
config_file = Path(__file__).parent / "configBoundaryLayer.cfg"
output_path = Path(__file__).parent / Path("output")
output_path.mkdir(exist_ok=True)


# %%
logger.Logger.initialize_logger()

# %%
communicator = parallel.MPICommunicator.get_instance()

config = basics.ConfigurationFile()
config.load(str(config_file))

para = gpu.Parameter(communicator.get_number_of_processes(), communicator.get_process_id(), config)
bc_factory = gpu.BoundaryConditionFactory()

grid_scaling_factory = gpu.GridScalingFactory()
grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)

# %%
boundary_layer_height = config.get_float_value("boundaryLayerHeight", 1000)
z0 = config.get_float_value("z0", 0.1)
u_star = config.get_float_value("u_star", 0.4)

von_karman_constant = config.get_float_value("vonKarmanConstant", 0.4)  # von Karman constant

viscosity = config.get_float_value("viscosity", 1.56e-5)

velocity = (
    0.5 * u_star / von_karman_constant * np.log(boundary_layer_height / z0)
)  # 0.5 times max mean velocity at the top in m/s

mach = config.get_float_value("Ma", 0.1)
nodes_per_height = config.get_uint_value("nz", 64)


write_precursor = config.get_bool_value("writePrecursor", False)
read_precursor = config.get_bool_value("readPrecursor", False)

# all in s
t_start_out = config.get_float_value("tStartOut")
t_out = config.get_float_value("tOut")
t_end = config.get_float_value("tEnd")  # total time of simulation

t_start_averaging = config.get_float_value("tStartAveraging")
t_start_tmp_averaging = config.get_float_value("tStartTmpAveraging")
t_averaging = config.get_float_value("tAveraging")
t_start_out_probe = config.get_float_value("tStartOutProbe")
t_out_probe = config.get_float_value("tOutProbe")

# %%
length = np.array([6, 4, 1]) * boundary_layer_height
dx = boundary_layer_height / nodes_per_height
dt = dx * mach / (np.sqrt(3) * velocity)
velocity_LB = velocity * dt / dx  # LB units
viscosity_LB = viscosity * dt / (dx * dx)  # LB units
pressure_gradient = u_star * u_star / boundary_layer_height
pressure_gradient_LB = pressure_gradient * (dt * dt) / dx

logger.vf_log_info(f"velocity  [dx/dt] = {velocity_LB}")
logger.vf_log_info(f"dt   = {dt}")
logger.vf_log_info(f"dx   = {dx}")
logger.vf_log_info(f"viscosity [10^8 dx^2/dt] = {viscosity_LB * 1e8}")
logger.vf_log_info(f"u* /(dx/dt) = {u_star * dt / dx}")
logger.vf_log_info(f"dpdx  = {pressure_gradient}")
logger.vf_log_info(f"dpdx /(dx/dt^2) = {pressure_gradient_LB}")

# %%
para.set_output_prefix(sim_name)
para.set_print_files(True)

para.set_forcing(pressure_gradient_LB, 0, 0)
para.set_velocity_LB(velocity_LB)
para.set_viscosity_LB(viscosity_LB)
para.set_velocity_ratio(dx / dt)
para.set_viscosity_ratio(dx * dx / dt)
para.set_density_ratio(1.0)

para.configure_main_kernel(gpu.kernel.compressible.K17CompressibleNavierStokes)

para.set_timestep_start_out(int(t_start_out / dt))
para.set_timestep_out(int(t_out / dt))
para.set_timestep_end(int(t_end / dt))
para.set_is_body_force(config.get_bool_value("bodyForce"))
# %%
tm_factory = gpu.TurbulenceModelFactory(para)
tm_factory.read_config_file(config)
# %%
grid_dimesions = gpu.grid_generator.GridDimensions(0, length[0], 0, length[1], 0, length[2], dx)
grid_builder = gpu.grid_generator.MultipleGridBuilderFacade(grid_dimesions)
grid_builder.set_periodic_boundary_condition(not read_precursor, True, False)
grid_builder.create_grids(communicator.get_process_id())

sampling_offset = 2
if read_precursor:
    precursor_directory = config.get_string_value("precursorDirectory")
    nTReadPrecursor = config.get_int_value("nTimestepsReadPrecursor")

    precursor = gpu.create_file_collection(precursor_directory, "precursor", gpu.TransientBCFileType.VTK)
    grid_builder.set_precursor_boundary_condition(gpu.SideType.MX, precursor, nTReadPrecursor, 0, 0, 0)

grid_builder.set_stress_boundary_condition(gpu.SideType.MZ, 0, 0, 1, sampling_offset, von_karman_constant, z0, dx)
grid_builder.set_slip_boundary_condition(gpu.SideType.PZ, 0, 0, -1)

if read_precursor:
    grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 0)
    bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)
    bc_factory.set_precursor_boundary_condition(gpu.PrecursorBC.PrecursorDistributions)

bc_factory.set_stress_boundary_condition(gpu.StressBC.StressBounceBackCompressible)
bc_factory.set_slip_boundary_condition(gpu.SlipBC.SlipCompressible)
para.set_outflow_pressure_correction_factor(0.0)
# %%
para.set_initial_condition_perturbed_log_law(u_star, z0, length[0], length[2], boundary_layer_height, dx / dt)
cuda_memory_manager = gpu.CudaMemoryManager(para)

# %%
planar_average_probe = gpu.probes.planar_average_probe.PlanarAverageProbe(
    para,
    cuda_memory_manager,
    para.get_output_path(),
    "horizontalPlanes",
    0,
    int(t_start_tmp_averaging / dt),
    int(t_averaging / dt),
    int(t_start_out_probe / dt),
    int(t_out_probe / dt),
    basics.geometry3d.Axis.z,
    True,
    False,
)
planar_average_probe.add_all_available_statistics()
planar_average_probe.set_file_name_to_n_out()
para.add_sampler(planar_average_probe)
# %%
wall_model_probe = gpu.probes.wall_model_probe.WallModelProbe(
    para,
    cuda_memory_manager,
    para.get_output_path(),
    "wallModelProbe",
    0,
    int(t_start_tmp_averaging / dt),
    int(t_averaging / dt / 4),
    int(t_start_out_probe / dt),
    int(t_out_probe / dt),
    False,
    True,
    True,
    para.get_is_body_force(),
    False,
)
para.add_sampler(wall_model_probe)

plane_locs = [
    100,
]
if read_precursor:
    plane_locs.extend([1000, 1500, 2000, 2500, 0])

for n_probe, probe_pos in enumerate(plane_locs):
    plane_probe = gpu.probes.probe.Probe(
        para,
        cuda_memory_manager,
        para.get_output_path(),
        f"planeProbe_{n_probe + 1}",
        int(t_start_averaging / dt),
        10,
        int(t_start_out_probe / dt),
        int(t_out_probe / dt),
        False,
        False,
    )
    plane_probe.set_probe_plane(probe_pos, 0, 0, dx, length[1], length[2])
    plane_probe.add_all_available_statistics()
    para.add_sampler(plane_probe)

if write_precursor:
    nTWritePrecursor = config.get_int_value("nTimestepsWritePrecursor")
    t_start_precursor = config.get_float_value("tStartPrecursor")
    pos_x_precursor = config.get_float_value("posXPrecursor")
    precursor_directory = config.get_string_value("precursorDirectory")

    precursor_writer = gpu.PrecursorWriter(
        para,
        cuda_memory_manager,
        para.get_output_path() + precursor_directory,
        "precursor",
        pos_x_precursor,
        0,
        length[1],
        0,
        length[2],
        int(t_start_precursor / dt),
        nTWritePrecursor,
        gpu.OutputVariable.Distributions,
        10000,
    )
    para.add_sampler(precursor_writer)
# %%
sim = gpu.Simulation(para, grid_builder.get_grid_builder(), bc_factory, tm_factory, grid_scaling_factory)
# %%
sim.init_timers()
for t in range(1, para.get_timestep_end() + 1):
    sim.calculate_timestep(t)
sim.finalize()
MPI.Finalize()
