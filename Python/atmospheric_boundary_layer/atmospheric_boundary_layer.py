#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger
#%%
boundary_layer_height = 2000.

length = np.array([6,4,1])*boundary_layer_height
viscosity = 1.56e-5
velocity = 9
mach = 0.1
nodes_per_height = 64

sim_name = "ABL"
config_file = Path(__file__).parent/Path("config.txt")
output_path = Path(__file__).parent/Path("output")
output_path.mkdir(exist_ok=True)
timeStepOut = 500
t_end = 1000

#%%
logger.Logger.initialize_logger()
basics.logger.Logger.add_stdout()
basics.logger.Logger.set_debug_level(basics.logger.Level.INFO_LOW)
basics.logger.Logger.time_stamp(basics.logger.TimeStamp.ENABLE)
basics.logger.Logger.enable_printed_rank_numbers(True)
#%%
grid_factory = gpu.grid.GridFactory.make()
grid_builder = gpu.grid.MultipleGridBuilder.make_shared(grid_factory)
dx = boundary_layer_height/nodes_per_height

grid_builder.add_coarse_grid(0.0, 0.0, 0.0, *length, dx)
grid_builder.set_periodic_boundary_condition(True, True, False)
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
para.set_main_kernel("TurbulentViscosityCumulantK17CompChim")
para.set_use_AMD(True)
para.set_SGS_constant(0.083)

def init_func(coord_x, coord_y, coord_z):
    return [0.0, velocity_lb, 0.0, 0.0]

para.set_initial_condition(init_func)
para.set_t_out(timeStepOut)
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
cuda_memory_manager = gpu.CudaMemoryManager.make(para)
grid_generator = gpu.grid.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager)
#%%
x_pos = length[0]/2
precursor_writer = gpu.PrecursorWriter("Precursor", str(output_path), x_pos, 0., length[1], 0., length[2], 50, 200, 5)
para.add_probe(precursor_writer)
#%%
#%%
sim = gpu.Simulation(comm)
kernel_factory = gpu.KernelFactory.get_instance()
sim.set_factories(kernel_factory, gpu.PreProcessorFactory.get_instance())
sim.init(para, grid_generator, gpu.FileWriter(), cuda_memory_manager)
#%%
sim.run()
sim.free()
MPI.Finalize()