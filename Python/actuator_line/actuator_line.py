#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger
#%%
reference_diameter = 126

length = np.array([30,8,8])*reference_diameter
viscosity = 1.56e-5
velocity = 9
mach = 0.1
nodes_per_diameter = 32

sim_name = "ActuatorLine"
config_file = Path(__file__).parent/Path("config.txt")
output_path = Path(__file__).parent/Path("output")
output_path.mkdir(exist_ok=True)
timeStepOut = 500
t_end = 50

#%%
logger.Logger.initialize_logger()
basics.logger.Logger.add_stdout()
basics.logger.Logger.set_debug_level(basics.logger.Level.INFO_LOW)
basics.logger.Logger.time_stamp(basics.logger.TimeStamp.ENABLE)
basics.logger.Logger.enable_printed_rank_numbers(True)
#%%
grid_builder = gpu.MultipleGridBuilder.make_shared()
dx = reference_diameter/nodes_per_diameter

grid_builder.add_coarse_grid(0.0, 0.0, 0.0, *length, dx)
grid_builder.set_periodic_boundary_condition(False, False, False)
grid_builder.build_grids(basics.LbmOrGks.LBM, False)
# %%
comm = gpu.Communicator.get_instance()
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
para.set_main_kernel("CumulantK17CompChim")

def init_func(coord_x, coord_y, coord_z):
    return [0.0, velocity_lb, 0.0, 0.0]

para.set_initial_condition(init_func)
para.set_t_out(timeStepOut)
para.set_t_end(int(t_end/dt))
para.set_is_body_force(True)

#%%
grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_lb, 0.0, 0.0)
grid_builder.set_velocity_boundary_condition(gpu.SideType.PX, velocity_lb, 0.0, 0.0)

grid_builder.set_velocity_boundary_condition(gpu.SideType.MY, velocity_lb, 0.0, 0.0)
grid_builder.set_velocity_boundary_condition(gpu.SideType.PY, velocity_lb, 0.0, 0.0)

grid_builder.set_velocity_boundary_condition(gpu.SideType.MZ, velocity_lb, 0.0, 0.0)
grid_builder.set_velocity_boundary_condition(gpu.SideType.PZ, velocity_lb, 0.0, 0.0)

#%%
cuda_memory_manager = gpu.CudaMemoryManager.make(para)
grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager)
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
point_probe = gpu.probes.PointProbe("pointProbe", str(output_path), 100, 500, 100)
point_probe.add_probe_points_from_list(np.array([1,2,5])*reference_diameter, np.array([3,3,3])*reference_diameter, np.array([3,3,3])*reference_diameter)
para.add_probe(point_probe)

plane_probe = gpu.probes.PlaneProbe("plane_probe", str(output_path), 100, 500, 100)
plane_probe.set_probe_plane(5*reference_diameter, 0, 0, dx, length[1], length[2])
para.add_probe(plane_probe)
#%%
sim = gpu.Simulation(comm)
kernel_factory = gpu.KernelFactory.get_instance()
sim.set_factories(kernel_factory, gpu.PreProcessorFactory.get_instance())
sim.init(para, grid_generator, gpu.FileWriter(), cuda_memory_manager)
#%%
sim.run()
sim.free()
MPI.Finalize()