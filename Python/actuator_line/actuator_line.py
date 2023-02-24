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

  You should have received a copy of the GNU General Public License along
  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.

! \file actuator_line.py
! \ingroup actuator_line
! \author Henry Korb, Henrik Asmuth
=======================================================================================
"""
#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids.bindings import basics, gpu, logger
#%%
sim_name = "ABL"
config_file = Path(__file__).parent/"configActuatorLine.txt"
output_path = Path(__file__).parent/Path("output")
output_path.mkdir(exist_ok=True)


#%%
logger.Logger.initialize_logger()
basics.logger.Logger.add_stdout()
basics.logger.Logger.set_debug_level(basics.logger.Level.INFO_LOW)
basics.logger.Logger.time_stamp(basics.logger.TimeStamp.ENABLE)
basics.logger.Logger.enable_printed_rank_numbers(True)
#%%
grid_factory = gpu.grid_generator.GridFactory.make()
grid_builder = gpu.grid_generator.MultipleGridBuilder.make_shared(grid_factory)
communicator = gpu.Communicator.get_instance()

config = basics.ConfigurationFile()
config.load(str(config_file))

para = gpu.Parameter(communicator.get_number_of_process(), communicator.get_pid(), config)
bc_factory = gpu.BoundaryConditionFactory()

#%%
turbine_diameter = config.get_float_value("turbineDiameter", 126)
boundary_layer_height = config.get_float_value("boundaryLayerHeight", 1000)
z0 = config.get_float_value("z0", 0.1)
u_star = config.get_float_value("u_star", 0.4)

kappa = config.get_float_value("vonKarmanConstant", 0.4) # von Karman constant

viscosity = config.get_float_value("viscosity", 1.56e-5)

velocity  = 0.5*u_star/kappa*np.log(boundary_layer_height/z0+1) #0.5 times max mean velocity at the top in m/s

mach = config.get_float_value("Ma", 0.1)
nodes_per_height = config.get_uint_value("nz", 64)


turb_pos = np.array([3,3,3])*turbine_diameter
epsilon = config.get_float_value("SmearingWidth", 5)
density = config.get_float_value("Density", 1.225)
level = 0
n_blades = 3
n_blade_nodes = config.get_int_value("NumberOfNodesPerAL", 32)

read_precursor = config.get_bool_value("readPrecursor", False)

if read_precursor:
    nTReadPrecursor = config.get_int_value("nTimestepsReadPrecursor")
    use_distributions = config.get_bool_value("useDistributions", False)
    precursor_directory = config.get_string_value("precursorDirectory")

# all in s
t_start_out   = config.get_float_value("tStartOut")
t_out        = config.get_float_value("tOut")
t_end        = config.get_float_value("tEnd") # total time of simulation

t_start_averaging     =  config.get_float_value("tStartAveraging")
t_start_tmp_averaging  =  config.get_float_value("tStartTmpAveraging")
t_averaging          =  config.get_float_value("tAveraging")
t_start_out_probe      =  config.get_float_value("tStartOutProbe")
t_out_probe           =  config.get_float_value("tOutProbe")

#%%
length = np.array([6,4,1])*boundary_layer_height
dx = boundary_layer_height/nodes_per_height
dt = dx * mach / (np.sqrt(3) * velocity)
velocity_ratio = dx/dt
velocity_LB = velocity / velocity_ratio  # LB units
viscosity_LB = viscosity / (velocity_ratio * dx)  # LB units
pressure_gradient = u_star * u_star / boundary_layer_height
pressure_gradient_LB = pressure_gradient * (dt*dt)/dx

logger.vf_log_info(f"velocity  [dx/dt] = {velocity_LB}")
logger.vf_log_info(f"dt   = {dt}")
logger.vf_log_info(f"dx   = {dx}")
logger.vf_log_info(f"viscosity [10^8 dx^2/dt] = {viscosity_LB*1e8}")
logger.vf_log_info(f"u* /(dx/dt) = {u_star*dt/dx}")
logger.vf_log_info(f"dpdx  = {pressure_gradient}")
logger.vf_log_info(f"dpdx /(dx/dt^2) = {pressure_gradient_LB}")


#%%
para.set_output_prefix(sim_name)
para.set_print_files(True)

para.set_forcing(pressure_gradient_LB, 0, 0)
para.set_velocity_LB(velocity_LB)
para.set_viscosity_LB(viscosity_LB)    
para.set_velocity_ratio(dx/dt)
para.set_viscosity_ratio(dx*dx/dt)
para.set_density_ratio(1.0)

para.set_main_kernel("CumulantK17")

para.set_timestep_start_out(int(t_start_out/dt))
para.set_timestep_out(int(t_out/dt))
para.set_timestep_end(int(t_end/dt))
para.set_is_body_force(True)
#%%
tm_factory = gpu.TurbulenceModelFactory(para)
tm_factory.read_config_file(config)
#%%
grid_scaling_factory = gpu.GridScalingFactory()
grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)

grid_builder.add_coarse_grid(0.0, 0.0, 0.0, length[0], length[1], length[2], dx)
grid_builder.set_periodic_boundary_condition(not read_precursor, True, False)
grid_builder.build_grids(basics.LbmOrGks.LBM, False)

sampling_offset = 2
if read_precursor:
    precursor = gpu.create_file_collection(precursor_directory + "/precursor", gpu.FileType.VTK)
    grid_builder.set_precursor_boundary_condition(gpu.SideType.MX, precursor, nTReadPrecursor, 0, 0, 0)

grid_builder.set_stress_boundary_condition(gpu.SideType.MZ, 0, 0, 1, sampling_offset, z0, dx)
para.set_has_wall_model_monitor(True)
grid_builder.set_slip_boundary_condition(gpu.SideType.PZ, 0, 0, -1)

if read_precursor:
    grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 0)
bc_factory.set_stress_boundary_condition(gpu.StressBC.StressPressureBounceBack)
bc_factory.set_slip_boundary_condition(gpu.SlipBC.SlipBounceBack) 
bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)
if read_precursor:
    bc_factory.set_precursor_boundary_condition(gpu.PrecursorBC.DistributionsPrecursor if use_distributions else gpu.PrecursorBC.VelocityPrecursor)
para.set_outflow_pressure_correction_factor(0.0)
#%%
# don't use python init functions, they are very slow! Just kept as an example.
# Define lambda in bindings and set it here.
# def init_func(coord_x, coord_y, coord_z):
#     return [
#         0.0, 
#         (u_star/0.4 * np.log(np.maximum(coord_z,z0)/z0) + 2.0*np.sin(np.pi*16*coord_x/length[0])*np.sin(np.pi*8*coord_z/boundary_layer_height)/(np.square(coord_z/boundary_layer_height)+1))  * dt / dx, 
#         2.0*np.sin(np.pi*16.*coord_x/length[0])*np.sin(np.pi*8.*coord_z/boundary_layer_height)/(np.square(coord_z/boundary_layer_height)+1.)  * dt / dx, 
#         8.0*u_star/0.4*(np.sin(np.pi*8.0*coord_y/boundary_layer_height)*np.sin(np.pi*8.0*coord_z/boundary_layer_height)+np.sin(np.pi*8.0*coord_x/length[0]))/(np.square(length[2]/2.0-coord_z)+1.) * dt / dx]
# para.set_initial_condition(init_func)
para.set_initial_condition_perturbed_log_law(u_star, z0, length[0], length[2], boundary_layer_height, velocity_ratio)

#%%
turb_pos = np.array([3, 3, 3])*turbine_diameter
epsilon = 1.5*dx
density = 1.225
level = 0
n_blades = 3
n_blade_nodes = 32
omega = 1
blade_radii = np.arange(n_blade_nodes, dtype=np.float32)/(0.5*turbine_diameter)
alm = gpu.ActuatorFarm(n_blades, density, n_blade_nodes, epsilon, level, dt, dx, True)
alm.add_turbine(turb_pos[0], turb_pos[1], turb_pos[2], turbine_diameter, omega, 0, 0, blade_radii)
para.add_actuator(alm)
#%%
planar_average_probe = gpu.probes.PlanarAverageProbe("horizontalPlanes", para.get_output_path(), 0, int(t_start_tmp_averaging/dt), int(t_averaging/dt) , int(t_start_out_probe/dt), int(t_out_probe/dt), 'z')
planar_average_probe.add_all_available_statistics()
planar_average_probe.set_file_name_to_n_out()
para.add_probe(planar_average_probe)
#%%
wall_model_probe = gpu.probes.WallModelProbe("wallModelProbe", para.get_output_path(), 0, int(t_start_tmp_averaging/dt), int(t_averaging/dt/4), int(t_start_out_probe/dt), int(t_out_probe/dt))
wall_model_probe.add_all_available_statistics()
wall_model_probe.set_file_name_to_n_out()
wall_model_probe.set_force_output_to_stress(True)
if para.get_is_body_force():
    wall_model_probe.set_evaluate_pressure_gradient(True)
para.add_probe(wall_model_probe)

plane_locs = [100,]
if read_precursor: plane_locs.extend([1000, 1500, 2000, 2500, 0])

for n_probe, probe_pos in enumerate(plane_locs):
    plane_probe = gpu.probes.PlaneProbe(f"planeProbe_{n_probe+1}", para.get_output_path(), int(t_start_averaging/dt), 10, int(t_start_out_probe/dt), int(t_out_probe/dt))
    plane_probe.set_probe_plane(probe_pos, 0, 0, dx, length[1], length[2])
    plane_probe.add_all_available_statistics()
    para.add_probe(plane_probe)
#%%
cuda_memory_manager = gpu.CudaMemoryManager(para)
grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager, communicator)
#%%
#%%
sim = gpu.Simulation(para, cuda_memory_manager, communicator, grid_generator, bc_factory, tm_factory, grid_scaling_factory)
#%%
sim.run()
MPI.Finalize()

