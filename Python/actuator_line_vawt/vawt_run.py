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
import numpy as np
import subprocess, sys
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger, parallel
from vawt_postprocess import load_wake, plot_wake, load_loads, plot_loads
from pyfluids.basics import ConfigurationFile
# define file paths
PATH_SCRIPT = Path(__file__).resolve().parent
PATH_DATA   = PATH_SCRIPT / "data"
PATH_RAW    = PATH_DATA / "output"
PATH_AIRFOIL = PATH_DATA / "airfoil"
PATH_POST   = PATH_DATA / "post"

#%% delete old data
subprocess.run(f"rm -rf {PATH_RAW}", shell=True, check=True)
subprocess.run(f"rm -rf {PATH_POST}", shell=True, check=True)
PATH_POST.mkdir(parents=True, exist_ok=True)
if len(sys.argv) < 2: raise ValueError("Please provide a config file, e.g. python3 vawt_run.py vawt_config_performance.cfg")
config_file = sys.argv[1]
logger.vf_log_info(f"using config: {config_file}")
PATH_CONFIG = PATH_SCRIPT / config_file
logger.Logger.initialize_logger() 

#%%
# initialize pyfluids modules
grid_builder    = gpu.grid_generator.MultipleGridBuilder()
communicator    = parallel.MPICommunicator.get_instance()
is_root_rank    = communicator.get_process_id() == 0
config          = basics.ConfigurationFile()
config.load(str(PATH_CONFIG))
para            = gpu.Parameter(communicator.get_number_of_processes(), communicator.get_process_id(), config)
bc_factory      = gpu.BoundaryConditionFactory()

# read parameters
viscosity = config.get_float_value("viscosity")
density = config.get_float_value("density")
name_airfoil_data = config.get_string_value("name_airfoil_data")
tip_speed_ratio = config.get_float_value("tip_speed_ratio")
number_blades = config.get_int_value("number_of_blades")
rotor_diameter = config.get_float_value("rotor_diameter")
rotor_height = config.get_float_value("rotor_height")
velocity_inlet = config.get_float_value("velocity_inlet")
reynolds_number = config.get_float_value("reynolds_number_chord")
blade_chord = config.get_float_value("blade_chord")
blade_pitch = config.get_float_value("blade_pitch_in_degrees")
blade_mounting_point = config.get_float_value("blade_mounting_position_chordwise")
number_cells_per_diameter_finest_level = config.get_int_value("number_of_cells_per_diameter")
number_of_points_per_blade = config.get_int_value("number_of_actuator_line_points_per_blade")
smearing_width = config.get_float_value("smearing_width_per_dx_fine")
flag_localized_smearing = config.get_bool_value("flag_localized_smearing_width")
number_fine_grid_levels = config.get_int_value("number_of_grid_refinement_levels")
quadric_delimiter = config.get_float_value("quadric_limiter")
position_domain_0_x_min = config.get_float_value("position_domain_0_x_min")
position_domain_0_y_min = config.get_float_value("position_domain_0_y_min")
position_domain_0_z_min = config.get_float_value("position_domain_0_z_min")
position_domain_0_x_max = config.get_float_value("position_domain_0_x_max")
position_domain_0_y_max = config.get_float_value("position_domain_0_y_max")
position_domain_0_z_max = config.get_float_value("position_domain_0_z_max")
position_domain_1_x_min = config.get_float_value("position_domain_1_x_min")
position_domain_1_y_min = config.get_float_value("position_domain_1_y_min")
position_domain_1_z_min = config.get_float_value("position_domain_1_z_min")
position_domain_1_x_max = config.get_float_value("position_domain_1_x_max")
position_domain_1_y_max = config.get_float_value("position_domain_1_y_max")
position_domain_1_z_max = config.get_float_value("position_domain_1_z_max")
position_turbine_x = config.get_float_value("position_turbine_x")
position_turbine_y = config.get_float_value("position_turbine_y")
position_turbine_z = config.get_float_value("position_turbine_z")
TurbulenceModel = config.get_string_value("TurbulenceModel")
SGSconstant = config.get_float_value("SGSconstant")
flag_flow_curvature = config.get_bool_value("flag_flow_curvature")
flag_end_effects = config.get_bool_value("flag_end_effects")
mach_number = config.get_float_value("mach_number_lbm")
flag_print_files = config.get_bool_value("flag_print_files")
nrev_end = config.get_float_value("number_of_revolutions_simulation_end")
position_probe_xplane_x_raw = config.get_string_value("position_probe_xplane_x")
position_probe_xplane_x = np.fromstring(position_probe_xplane_x_raw,sep=",",dtype=float)
position_probe_xplane_y_min = config.get_float_value("position_probe_xplane_y_min")
position_probe_xplane_z_min = config.get_float_value("position_probe_xplane_z_min")
position_probe_xplane_y_max = config.get_float_value("position_probe_xplane_y_max")
position_probe_xplane_z_max = config.get_float_value("position_probe_xplane_z_max")
position_probe_zplane_x_min = config.get_float_value("position_probe_zplane_x_min")
position_probe_zplane_y_min = config.get_float_value("position_probe_zplane_y_min")
position_probe_zplane_x_max = config.get_float_value("position_probe_zplane_x_max")
position_probe_zplane_y_max = config.get_float_value("position_probe_zplane_y_max")
position_probe_zplane_z = config.get_float_value("position_probe_zplane_z")
def _get_float_with_auto_fallback(config: ConfigurationFile, float_target_name: str, float_fallback: float) -> float:
    """
    provides float from config and directs provided
    fallback value as output in case the targeted
    parameter is defined as "auto" in the config file.
    """
    try: 
        return config.get_float_value(float_target_name)
    except RuntimeError:
        val = config.get_string_value(float_target_name)
        if val.lower() == "auto":
            return float_fallback
        else:
            raise ValueError(f"Invalid value for {float_target_name}: {val}")
nrev_start_write_loads = _get_float_with_auto_fallback(config,"number_of_revolutions_loads_output_start",nrev_end)
nrev_intervall_write_load = config.get_float_value("number_of_revolutions_loads_output_intervall")
nrev_start_write_probe = _get_float_with_auto_fallback(config,"number_of_revolutions_probe_output_start",nrev_end)
nrev_intervall_write_probe = config.get_float_value("number_of_revolutions_probe_output_intervall")
nrev_start_probe_averaging = _get_float_with_auto_fallback(config,"number_of_revolutions_probe_averaging_start",nrev_end)
nrev_intervall_averaging = config.get_float_value("number_of_revolutions_probe_averaging_intervall")
nrev_start_write_field = _get_float_with_auto_fallback(config,"number_of_revolutions_field_output_start",nrev_end)
nrev_intervall_write_field = config.get_float_value("number_of_revolutions_field_output_intervall")

# vawt
rotor_speed = tip_speed_ratio * velocity_inlet / (rotor_diameter/2)  # Rotor angular velocity
time_revolution = 2 * np.pi / rotor_speed  # Time for one full revolution

# lbm
number_cells_per_diameter_coarse_level = number_cells_per_diameter_finest_level / (2**number_fine_grid_levels)
dx = rotor_diameter / number_cells_per_diameter_coarse_level     
dt = mach_number * (dx/(np.sqrt(3)*velocity_inlet))
velocity_ratio = dx/dt  
velocity_lb = velocity_inlet / velocity_ratio
viscosity_lb = viscosity / (velocity_ratio * dx)  

# airfoil
airfoil_data = np.loadtxt(PATH_AIRFOIL / f"{name_airfoil_data}.dat", comments="%", usecols=(0, 1, 2))
polar_alpha = airfoil_data[:,0]
polar_cl    = airfoil_data[:,1]
polar_cd    = airfoil_data[:,2]

# compute physical sampling times from revolution specific sampling times (as set in configuration file)
time_start_write_loads      = time_revolution * nrev_start_write_loads
time_start_write_probe      = time_revolution * nrev_start_write_probe
time_start_write_field      = time_revolution * nrev_start_write_field
time_start_probe_averaging  = time_revolution * nrev_start_probe_averaging
time_intervall_write_loads  = max(time_revolution * nrev_intervall_write_load, dt)
time_intervall_write_probe  = max(time_revolution * nrev_intervall_write_probe, dt)
time_intervall_write_field  = max(time_revolution * nrev_intervall_write_field, dt)
time_intervall_averaging    = time_revolution * nrev_intervall_averaging
time_end                    = time_revolution * nrev_end
number_steps_per_revolution = time_revolution / dt

# setup
para.set_output_path(str(PATH_RAW))
para.set_output_prefix("fields/vawt")
para.set_print_files(flag_print_files)
para.set_velocity_LB(velocity_lb)
para.set_viscosity_LB(viscosity_lb)
para.set_velocity_ratio(dx / dt)
para.set_viscosity_ratio(dx * dx / dt)
para.set_density_ratio(1.0)
para.configure_main_kernel(gpu.kernel.compressible.K17CompressibleNavierStokes)
para.set_timestep_start_out(int(time_start_write_field / dt))
para.set_timestep_out(int(time_intervall_write_field/dt))
para.set_timestep_end(int(time_end / dt))
para.set_is_body_force(True)
para.set_quadric_limiters(quadric_delimiter,quadric_delimiter,quadric_delimiter)

# logging
logger.vf_log_info(f"time_start_write_probe / dt = {int(time_start_write_probe / dt )}")
logger.vf_log_info(f"time_start_probe_averaging / dt = {int(time_start_probe_averaging / dt )}")
logger.vf_log_info(f"time_start_write_field / dt = {int(time_start_write_field / dt )}")
logger.vf_log_info(f"time_intervall_write_probe / dt = {int(time_intervall_write_probe / dt )}")
logger.vf_log_info(f"time_intervall_averaging / dt = {int(time_intervall_averaging / dt )}")
logger.vf_log_info(f"time_intervall_write_field / dt = {int(time_intervall_write_field / dt )}")
logger.vf_log_info(f"time_end / dt = {int(time_end / dt )}")
logger.vf_log_info(f"number_steps_per_revolution = {number_steps_per_revolution}")
logger.vf_log_info(f"time for one revolution = {time_revolution:.4f} s")
logger.vf_log_info(f"lbm velocity [dx/dt] = {velocity_lb}")
logger.vf_log_info(f"dt = {dt}")
logger.vf_log_info(f"dx = {dx}")
logger.vf_log_info(f"lbm viscosity [10^8 dx^2/dt] = {viscosity_lb * 1e8}")
logger.vf_log_info(f"rotor_speed [rad/s] = {rotor_speed}")
logger.vf_log_info(f"mesh cells per turbine rotor_diameter on finest grid_level_for_alm = {number_cells_per_diameter_finest_level}")
logger.vf_log_info(f"mesh cells per turbine rotor_diameter on coarsest grid_level_for_alm = {number_cells_per_diameter_coarse_level}")
logger.vf_log_info(f"mach number  = {mach_number}")
logger.vf_log_info(f"vawt sanityt check: mach_number {mach_number} < {np.sqrt(3)*velocity_inlet/(tip_speed_ratio*velocity_inlet)}") 

# set turbulence model
tm_factory = gpu.TurbulenceModelFactory(para)
tm_factory.read_config_file(config)

# create mesh
grid_scaling_factory = gpu.GridScalingFactory()
grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)            
grid_builder.add_coarse_grid(position_domain_0_x_min,
                             position_domain_0_y_min,
                             position_domain_0_z_min,
                             position_domain_0_x_max,
                             position_domain_0_y_max,
                             position_domain_0_z_max,
                             dx)
dx_fine = dx
grid_level_for_alm = 0
if number_fine_grid_levels > 0 :
    grid_builder.add_grid(gpu.grid_generator.Cuboid(position_domain_1_x_min,
                                                    position_domain_1_y_min,
                                                    position_domain_1_z_min,
                                                    position_domain_1_x_max,
                                                    position_domain_1_y_max,
                                                    position_domain_1_z_max),number_fine_grid_levels)
    grid_level_for_alm = number_fine_grid_levels
    dx_fine = dx / (2 ** grid_level_for_alm)
grid_builder.build_grids(False)

# assign boundary conditions to each wall
grid_builder.set_periodic_boundary_condition(False, False, False)
grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_lb, 0, 0)
grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 1.0)
grid_builder.set_slip_boundary_condition(gpu.SideType.MZ, 0, 0, 1)
grid_builder.set_slip_boundary_condition(gpu.SideType.PZ, 0, 0, -1)
grid_builder.set_slip_boundary_condition(gpu.SideType.MY, 0, 1, 0)
grid_builder.set_slip_boundary_condition(gpu.SideType.PY, 0, -1, 0)

# define boundary condition variant for each type
bc_factory.set_velocity_boundary_condition(gpu.VelocityBC.VelocityWithPressureInterpolatedCompressible) # OPTIONS: VelocityInterpolatedCompressible, VelocityBounceBack, VelocityInterpolatedIncompressible, VelocityWithPressureInterpolatedCompressible
bc_factory.set_slip_boundary_condition(gpu.SlipBC.SlipTurbulentViscosityCompressible)
bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)
para.set_outflow_pressure_correction_factor(0.0)

# initial conditions
para.set_initial_condition_uniform(velocity_lb, 0.0, 0.0)

# load parameters into cuda memory manager
cuda_memory_manager = gpu.CudaMemoryManager(para)

# set up fail-back for smearing width
# smearing_width = max(smearing_width*dx_fine,blade_chord) # max(..,blade_chord) avoids numerical artefacts for fine resolutions (empirical fallback)
smearing_width = smearing_width * dx_fine

# set up actuator farm kernel
alm = gpu.ActuatorFarmStandaloneVAWT(
    para,
    cuda_memory_manager,
    rotor_diameter,
    number_blades,
    number_of_points_per_blade,
    rotor_height,
    position_turbine_x*np.ones(1), # np.ones(1) required, as the actuator farm class is defined for possible multiple turbines
    position_turbine_y*np.ones(1),
    position_turbine_z*np.ones(1),
    rotor_speed*np.ones(1), 
    smearing_width,
    grid_level_for_alm,
    polar_alpha,
    polar_cl,
    polar_cd,
    blade_chord,
    blade_pitch,
    blade_mounting_point,
    velocity_inlet,
    flag_localized_smearing,
    flag_flow_curvature,
    flag_end_effects,
)
alm.enable_output("loads/ALM", int(time_start_write_loads / dt), int(time_intervall_write_loads / dt))
para.add_interactor(alm)

# find probe plane positions that are closest to target plane position but also align with the next fine grid node plane. 
position_probe_xplane_x = np.round((position_probe_xplane_x - position_domain_0_x_min) / dx_fine) * dx_fine + position_domain_0_x_min
position_probe_zplane_z = round((position_probe_zplane_z - position_domain_0_z_min) / dx_fine) * dx_fine + position_domain_0_z_min

# add x-plane probes
for index, x_position in enumerate(position_probe_xplane_x):
    plane_probe_x = gpu.probes.probe.Probe(
        para,
        cuda_memory_manager,
        para.get_output_path(),
        str(Path("./planes") / f"x_plane_{index}"),
        max(int(time_start_probe_averaging / dt), 1),
        max(int(time_intervall_averaging / dt), 1),
        max(int(time_start_write_probe / dt), 1),
        max(int(time_intervall_write_probe / dt), 1),
        False, # output timeseries
        False, # average every timestep
    )
    # define probe geometry
    plane_probe_x.set_probe_plane(x_position,
                                position_probe_xplane_y_min,
                                position_probe_xplane_z_min,
                                dx_fine,
                                position_probe_xplane_y_max-position_probe_xplane_y_min,
                                position_probe_xplane_z_max-position_probe_xplane_z_min)
    # add staistics + add to sampler
    plane_probe_x.add_all_available_statistics()
    para.add_sampler(plane_probe_x)
# add z-plane probe
plane_probe_z = gpu.probes.probe.Probe(
    para,
    cuda_memory_manager,
    para.get_output_path(),
    str(Path("./planes") / "z_plane"),
    max(int(time_start_probe_averaging / dt), 1),
    max(int(time_intervall_averaging / dt), 1),
    max(int(time_start_write_probe / dt), 1),
    max(int(time_intervall_write_probe / dt), 1),
    False, # output timeseries
    False, # average every timestep
)
# define probe geometry
plane_probe_z.set_probe_plane(position_probe_zplane_x_min,
                            position_probe_zplane_y_min,
                            position_probe_zplane_z,
                            position_probe_zplane_x_max-position_probe_zplane_x_min,
                            position_probe_zplane_y_max-position_probe_zplane_y_min,
                            dx_fine)
# add staistics + add to sampler
plane_probe_z.add_all_available_statistics()
para.add_sampler(plane_probe_z)

# create simulation object
sim = gpu.Simulation(para, cuda_memory_manager, grid_builder, bc_factory, tm_factory, grid_scaling_factory)

#%%
sim.run()
MPI.Finalize()

#%% 
logger.vf_log_info(".... creating loads plots for Figure 5/7 in TORQUE paper")
alpha_dict, fn_dict, ft_dict = load_loads(number_blades, flag_flow_curvature, flag_end_effects)
plot_loads(alpha_dict, fn_dict, ft_dict, flag_flow_curvature, flag_end_effects, PATH_POST / "loads.png")
if flag_end_effects and flag_flow_curvature: # wake plots only relevant for validation case.
    logger.vf_log_info(".... creating wake plots for Figure 8 in TORQUE paper")
    vx_y_dict, vx_z_dict, vy_y_dict = load_wake(dx_fine)
    plot_wake(vx_y_dict, vx_z_dict, vy_y_dict, PATH_POST / "wake.png")