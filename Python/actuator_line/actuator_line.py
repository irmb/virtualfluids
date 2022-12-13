#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger
#%%
reference_diameter = 126

length = np.array([29,6,6])*reference_diameter
viscosity = 1.56e-5
velocity = 9
mach = 0.1
nodes_per_diameter = 32

sim_name = "ActuatorLine"
config_file = Path(__file__).parent/Path("config.txt")
output_path = Path(__file__).parent/Path("output")
output_path.mkdir(exist_ok=True)
t_out = 100.
t_end = 500.

#%%
logger.Logger.initialize_logger()
basics.logger.Logger.add_stdout()
basics.logger.Logger.set_debug_level(basics.logger.Level.INFO_LOW)
basics.logger.Logger.time_stamp(basics.logger.TimeStamp.ENABLE)
basics.logger.Logger.enable_printed_rank_numbers(True)
# %%
comm = gpu.Communicator.get_instance()
#%%
grid_factory = gpu.grid_generator.GridFactory.make()
grid_builder = gpu.grid_generator.MultipleGridBuilder.make_shared(grid_factory)

#%%
dx = reference_diameter/nodes_per_diameter

grid_builder.add_coarse_grid(0.0, 0.0, 0.0, *length, dx)
grid_builder.set_periodic_boundary_condition(False, False, False)
grid_builder.build_grids(basics.LbmOrGks.LBM, False)
#%%
config = basics.ConfigurationFile()
config.load(str(config_file))
#%%
para = gpu.Parameter(config, comm.get_number_of_process(), comm.get_pid())

dt = dx * mach / (np.sqrt(3) * velocity)
velocity_lb = velocity * dt / dx # LB units
viscosity_lb = viscosity * dt / (dx * dx) # LB units

#%%
para.set_devices([0])
para.set_output_prefix(sim_name)
para.set_output_path(str(output_path))
para.set_f_name(para.get_output_path() + "/" + para.get_output_prefix())
para.set_print_files(True)
para.set_max_level(1)
#%%
para.set_velocity(velocity_lb)
para.set_viscosity(viscosity_lb)    
para.set_velocity_ratio(dx/dt)
para.set_viscosity_ratio(dx*dx/dt)
para.set_main_kernel("TurbulentViscosityCumulantK17CompChim")
para.set_use_AMD(True)
para.set_SGS_constant(0.083)

def init_func(coord_x, coord_y, coord_z):
    return [0.0, velocity_lb, 0.0, 0.0]

para.set_initial_condition(init_func)
para.set_t_out(int(t_out/dt))
para.set_t_end(int(t_end/dt))
para.set_is_body_force(True)

#%%
grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_lb, 0.0, 0.0)

grid_builder.set_velocity_boundary_condition(gpu.SideType.MY, velocity_lb, 0.0, 0.0)
grid_builder.set_velocity_boundary_condition(gpu.SideType.PY, velocity_lb, 0.0, 0.0)

grid_builder.set_velocity_boundary_condition(gpu.SideType.MZ, velocity_lb, 0.0, 0.0)
grid_builder.set_velocity_boundary_condition(gpu.SideType.PZ, velocity_lb, 0.0, 0.0)

grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 0.0)

#%%
cuda_memory_manager = gpu.CudaMemoryManager(para)
grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager, comm)
#%%
turb_pos = np.array([3,3,3])*reference_diameter
epsilon = 5
density = 1.225
level = 0
n_blades = 3
n_blade_nodes = 32
alm = gpu.ActuatorLine(n_blades, density, n_blade_nodes, epsilon, *turb_pos, reference_diameter, level, dt, dx)
para.add_actuator(alm)
#%%
point_probe = gpu.probes.PointProbe("pointProbe", str(output_path), 100, 1, 500, 100)
point_probe.add_probe_points_from_list(np.array([1,2,5])*reference_diameter, np.array([3,3,3])*reference_diameter, np.array([3,3,3])*reference_diameter)
point_probe.add_statistic(gpu.probes.Statistic.Means)

para.add_probe(point_probe)

plane_probe = gpu.probes.PlaneProbe("planeProbe", str(output_path), 100, 1, 500, 100)
plane_probe.set_probe_plane(5*reference_diameter, 0, 0, dx, length[1], length[2])
para.add_probe(plane_probe)
#%%
sim = gpu.Simulation(para, cuda_memory_manager, comm, grid_generator)
#%%
sim.run()
MPI.Finalize()