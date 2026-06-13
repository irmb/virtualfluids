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
#%%
import re
from pathlib import Path
import numpy as np
from mpi4py import MPI
from pyfluids import basics, gpu, logger, parallel

LOAD_FILE_PATTERN = re.compile(r"^ALM_ID_(?P<process_id>\d+)_t_(?P<t>\d+)\.bin\.vtu$")
CENTERLINE_PROBE_PATTERN = re.compile(r"^xLine(?P<probe_id>\d+)$")

PATH_SCRIPT = Path(__file__).resolve().parent
PATH_DATA = PATH_SCRIPT / "data"
PATH_RAW = PATH_DATA / "output"
PATH_POST = PATH_DATA / "post"
PATH_CONFIG = PATH_SCRIPT / "disk_config.cfg"

#%% delete old data
import subprocess
subprocess.run(f"rm -rf {PATH_RAW}", shell=True, check=True)
subprocess.run(f"rm -rf {PATH_POST}", shell=True, check=True)
PATH_RAW.mkdir(parents=True, exist_ok=True)
PATH_POST.mkdir(parents=True, exist_ok=True)

#%%
communicator = parallel.MPICommunicator.get_instance()
config = basics.ConfigurationFile()
config.load(str(PATH_CONFIG))

# parameters from config 
viscosity = config.get_float_value("viscosity")
number_of_blades = config.get_int_value("number_of_blades")
disk_diameter = config.get_float_value("disk_diameter")
velocity_inlet = config.get_float_value("velocity_inlet")
number_of_thrust_coefficients = config.get_int_value("number_of_thrust_coefficients")
thrust_coefficients = np.asarray(
    [config.get_float_value(f"thrust_coefficient_{index}") for index in range(number_of_thrust_coefficients)],
    dtype=float,
)
number_of_cells_per_diameter = config.get_int_value("number_of_cells_per_diameter")
number_of_actuator_line_points_per_blade = config.get_int_value("number_of_actuator_line_points_per_blade")
smearing_width_per_dx = config.get_float_value("smearing_width_per_dx")
position_domain_0_x_min = config.get_float_value("position_domain_0_x_min")
position_domain_0_y_min = config.get_float_value("position_domain_0_y_min")
position_domain_0_z_min = config.get_float_value("position_domain_0_z_min")
position_domain_0_x_max = config.get_float_value("position_domain_0_x_max")
position_domain_0_y_max = config.get_float_value("position_domain_0_y_max")
position_domain_0_z_max = config.get_float_value("position_domain_0_z_max")
position_disk_center_x = config.get_float_value("position_disk_center_x")
position_disk_center_y = config.get_float_value("position_disk_center_y")
position_disk_center_z = config.get_float_value("position_disk_center_z")
quadric_limiter = config.get_float_value("quadric_limiter")
mach_number = config.get_float_value("mach_number_lbm")
flag_print_files = config.get_bool_value("flag_print_files")
number_of_flow_throughs_loads_output_start = config.get_float_value("number_of_flow_throughs_loads_output_start")
number_of_flow_throughs_probe_output_start = config.get_float_value("number_of_flow_throughs_probe_output_start")
number_of_flow_throughs_field_output_start = config.get_float_value("number_of_flow_throughs_field_output_start")
number_of_flow_throughs_probe_averaging_start = config.get_float_value("number_of_flow_throughs_probe_averaging_start")
number_of_flow_throughs_probe_averaging_intervall = config.get_float_value("number_of_flow_throughs_probe_averaging_intervall")
number_of_flow_throughs_loads_output_intervall = config.get_float_value("number_of_flow_throughs_loads_output_intervall")
number_of_flow_throughs_probe_output_intervall = config.get_float_value("number_of_flow_throughs_probe_output_intervall")
number_of_flow_throughs_field_output_intervall = config.get_float_value("number_of_flow_throughs_field_output_intervall")
number_of_flow_throughs_simulation_end = config.get_float_value("number_of_flow_throughs_simulation_end")
position_probe_zplane_x_min = config.get_float_value("position_probe_zplane_x_min")
position_probe_zplane_y_min = config.get_float_value("position_probe_zplane_y_min")
position_probe_zplane_x_max = config.get_float_value("position_probe_zplane_x_max")
position_probe_zplane_y_max = config.get_float_value("position_probe_zplane_y_max")
position_probe_zplane_z = config.get_float_value("position_probe_zplane_z")

# pre-process parameters
xcm = position_domain_0_x_min
ycm = position_domain_0_y_min
zcm = position_domain_0_z_min
xcp = position_domain_0_x_max
ycp = position_domain_0_y_max
zcp = position_domain_0_z_max
disk_center = np.asarray([position_disk_center_x, position_disk_center_y, position_disk_center_z], dtype=float)
position_disk = disk_center[:, np.newaxis]
time_flow_through_domain = (xcp - xcm) / velocity_inlet
dx = disk_diameter / number_of_cells_per_diameter
dt = mach_number * (dx / (np.sqrt(3.0) * velocity_inlet))
velocity_ratio = dx / dt
velocity_lb = velocity_inlet / velocity_ratio
viscosity_lb = viscosity / (velocity_ratio * dx)
smearing_width = smearing_width_per_dx * dx

if thrust_coefficients.size == 0:
    raise ValueError("At least one thrust coefficient must be defined.")

timestep_start_write_loads = max(int(time_flow_through_domain * number_of_flow_throughs_loads_output_start / dt), 1)
timestep_start_write_probe = max(int(time_flow_through_domain * number_of_flow_throughs_probe_output_start / dt), 1)
timestep_start_write_field = max(int(time_flow_through_domain * number_of_flow_throughs_field_output_start / dt), 1)
timestep_start_probe_averaging = max(int(time_flow_through_domain * number_of_flow_throughs_probe_averaging_start / dt), 1)
timestep_interval_averaging = max(int(time_flow_through_domain * number_of_flow_throughs_probe_averaging_intervall / dt), 1)
timestep_interval_write_loads = max(int(time_flow_through_domain * number_of_flow_throughs_loads_output_intervall / dt), 1)
timestep_interval_write_probe = max(int(time_flow_through_domain * number_of_flow_throughs_probe_output_intervall / dt), 1)
timestep_interval_write_field = max(int(time_flow_through_domain * number_of_flow_throughs_field_output_intervall / dt), 1)
timestep_end = int(time_flow_through_domain * number_of_flow_throughs_simulation_end / dt)

logger.Logger.initialize_logger()
logger.vf_log_info(f"time_start_write_probe / dt = {timestep_start_write_probe}")
logger.vf_log_info(f"time_start_probe_averaging / dt = {timestep_start_probe_averaging}")
logger.vf_log_info(f"time_start_write_field / dt = {timestep_start_write_field}")
logger.vf_log_info(f"time_interval_write_probe / dt = {timestep_interval_write_probe}")
logger.vf_log_info(f"time_interval_averaging / dt = {timestep_interval_averaging}")
logger.vf_log_info(f"time_interval_write_field / dt = {timestep_interval_write_field}")
logger.vf_log_info(f"time_end / dt = {timestep_end}")
logger.vf_log_info(f"time for one flow-through = {time_flow_through_domain:.4f} s")
logger.vf_log_info(f"lbm velocity [dx/dt] = {velocity_lb}")
logger.vf_log_info(f"dt = {dt}")
logger.vf_log_info(f"dx = {dx}")
logger.vf_log_info(f"smearing width [m] = {smearing_width}")

#%% run for different thrust coefficients
for thrust_coefficient in thrust_coefficients :
    case_tag = f"ct_{thrust_coefficient:.2f}".replace(".", "p")
    case_raw_path = PATH_RAW / case_tag

    para = gpu.Parameter(communicator.get_number_of_processes(), communicator.get_process_id(), config)
    bc_factory = gpu.BoundaryConditionFactory()
    tm_factory = gpu.TurbulenceModelFactory(para)
    tm_factory.read_config_file(config)

    grid_builder = gpu.grid_generator.MultipleGridBuilder()
    grid_scaling_factory = gpu.GridScalingFactory()
    grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)

    grid_builder.add_coarse_grid(xcm, ycm, zcm, xcp, ycp, zcp, dx)
    grid_level_for_alm = 0
    grid_builder.build_grids(False)

    grid_builder.set_periodic_boundary_condition(False, False, False)
    grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_lb, 0.0, 0.0)
    grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 1.0)
    grid_builder.set_slip_boundary_condition(gpu.SideType.MZ, 0.0, 0.0, 1.0)
    grid_builder.set_slip_boundary_condition(gpu.SideType.PZ, 0.0, 0.0, -1.0)
    grid_builder.set_slip_boundary_condition(gpu.SideType.MY, 0.0, 1.0, 0.0)
    grid_builder.set_slip_boundary_condition(gpu.SideType.PY, 0.0, -1.0, 0.0)

    bc_factory.set_velocity_boundary_condition(gpu.VelocityBC.VelocityWithPressureInterpolatedCompressible)
    bc_factory.set_slip_boundary_condition(gpu.SlipBC.SlipTurbulentViscosityCompressible)
    bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)

    para.set_outflow_pressure_correction_factor(0.0)
    para.set_initial_condition_uniform(velocity_lb, 0.0, 0.0)
    para.set_output_path(str(case_raw_path))
    para.set_output_prefix("fields/disk")
    para.set_print_files(flag_print_files)
    para.set_velocity_LB(velocity_lb)
    para.set_viscosity_LB(viscosity_lb)
    para.set_velocity_ratio(velocity_ratio)
    para.set_viscosity_ratio(dx * dx / dt)
    para.set_density_ratio(1.0)
    para.configure_main_kernel(gpu.kernel.compressible.K17CompressibleNavierStokes)
    para.set_timestep_start_out(timestep_start_write_field)
    para.set_timestep_out(timestep_interval_write_field)
    para.set_timestep_end(timestep_end)
    para.set_is_body_force(True)
    para.set_quadric_limiters(quadric_limiter, quadric_limiter, quadric_limiter)

    cuda_memory_manager = gpu.CudaMemoryManager(para)

    # actuator disk quantities
    center_radius = 0.25 * disk_diameter / number_of_actuator_line_points_per_blade
    induction = 0.5 * (1.0 - np.sqrt(1.0 - float(thrust_coefficient))) # 1D momentum theory
    disk_normal_coefficient = 4.0 * induction / (1.0 - induction)
    hub_drag_coefficient = disk_normal_coefficient

    hub_config = gpu.HubConfig(0.0, center_radius, 0.0, 1)
    hub_skin_friction_coefficient = 0.0

    blade_radii = (0.5 * disk_diameter / number_of_actuator_line_points_per_blade) * (np.arange(number_of_actuator_line_points_per_blade, dtype=float) + 0.5)

    inner_radii = np.empty_like(blade_radii)
    outer_radii = np.empty_like(blade_radii)
    inner_radii[0] = center_radius
    inner_radii[1:] = 0.5 * (blade_radii[1:] + blade_radii[:-1])
    outer_radii[:-1] = 0.5 * (blade_radii[1:] + blade_radii[:-1])
    outer_radii[-1] = 0.5 * disk_diameter
    annulus_areas = np.pi * (outer_radii**2 - inner_radii**2) / number_of_blades

    previous_radii = np.zeros_like(blade_radii)
    previous_radii[1:] = blade_radii[:-1]
    next_radii = np.empty_like(blade_radii)
    next_radii[:-1] = blade_radii[1:]
    next_radii[-1] = 0.5 * disk_diameter
    blade_point_widths = 0.5 * (next_radii - previous_radii)

    blade_chords = 2.0 * np.sqrt(np.clip(1.0 - (4.0 * blade_radii / disk_diameter - 1.0) ** 2, 0.0, None))
    blade_normal_coefficients = (disk_normal_coefficient * annulus_areas / (blade_chords * blade_point_widths)).tolist()

    # configure
    alm = gpu.ActuatorFarmStandalone(
        para,
        cuda_memory_manager,
        disk_diameter,
        number_of_actuator_line_points_per_blade,
        [disk_center[0]],
        [disk_center[1]],
        [disk_center[2]],
        [0.0],
        smearing_width,
        grid_level_for_alm,
        hub_config=hub_config,
        hub_drag_coeff=hub_drag_coefficient,
        hub_skin_friction_coeff=hub_skin_friction_coefficient,
        number_of_blades=number_of_blades,
        blade_normal_coefficients=blade_normal_coefficients,
    )
    alm.enable_output("loads/ALM", timestep_start_write_loads, timestep_interval_write_loads)
    para.add_interactor(alm)

    # add probe plane
    x0, x1 = sorted((position_probe_zplane_x_min, position_probe_zplane_x_max))
    y0, y1 = sorted((position_probe_zplane_y_min, position_probe_zplane_y_max))
    z_position_requested = position_probe_zplane_z
    level0_min_z = zcm - 0.5 * dx
    lower = level0_min_z + np.floor((z_position_requested - level0_min_z) / dx) * dx
    upper = lower + dx
    z0 = float(lower if (z_position_requested - lower) <= (upper - z_position_requested) else upper)

    probe = gpu.probes.probe.Probe(
        para,
        cuda_memory_manager,
        para.get_output_path(),
        "zplane/zPlane",
        timestep_start_probe_averaging,
        timestep_interval_averaging,
        timestep_start_write_probe,
        timestep_interval_write_probe,
        False,
        False,
    )
    probe.set_probe_plane(x0, y0, z0, x1 - x0, y1 - y0, float(dx))
    probe.add_all_available_statistics()
    para.add_sampler(probe)

    # run simulation
    sim = gpu.Simulation(para, cuda_memory_manager, grid_builder, bc_factory, tm_factory, grid_scaling_factory)
    sim.run()

#%%
if communicator.get_process_id() == 0:
    from disk_postprocess import load_velocities, plot_velocities
    dict_vx = load_velocities(PATH_RAW, thrust_coefficients, velocity_inlet, position_disk, disk_diameter)
    plot_velocities(dict_vx, PATH_POST / "velocities.png", position_disk)
MPI.Finalize()
