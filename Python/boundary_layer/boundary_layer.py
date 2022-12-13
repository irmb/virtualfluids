#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger
#%%
reference_height = 1000 # boundary layer height in m

length = np.array([6,4,1])*reference_height
viscosity = 1.56e-5
mach = 0.1
nodes_per_height = 32

z_0 = 0.1
u_star = 0.4
kappa = 0.4

velocity = 0.5*u_star/kappa*np.log(length[2]/z_0+1)
flow_through_time = length[0]/velocity
use_AMD = True


sim_name = "BoundaryLayer"
config_file = Path(__file__).parent/Path("config.txt")
output_path = Path(__file__).parent/Path("output")
output_path.mkdir(exist_ok=True)
t_out = 1000.
t_end = 5000.

t_start_averaging = 0
t_start_tmp_averaging =  100_000
t_averaging = 200
t_start_out_probe = 0
t_out_probe = 1000

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
dx = reference_height/nodes_per_height
dt = dx * mach / (np.sqrt(3) * velocity)
velocity_lb = velocity * dt / dx # LB units
viscosity_lb = viscosity * dt / (dx * dx) # LB units

pressure_gradient = u_star**2 / reference_height
pressure_gradient_lb = pressure_gradient * dt**2 / dx

logger.vf_log_info(f"velocity    = {velocity_lb:1.6} dx/dt")
logger.vf_log_info(f"dt          = {dt:1.6}")
logger.vf_log_info(f"dx          = {dx:1.6}")
logger.vf_log_info(f"u*          = {u_star:1.6}")
logger.vf_log_info(f"dpdx        = {pressure_gradient:1.6}")
logger.vf_log_info(f"dpdx        = {pressure_gradient_lb:1.6} dx/dt^2")
logger.vf_log_info(f"viscosity   = {viscosity_lb:1.6} dx^2/dt")


#%%
config = basics.ConfigurationFile()
config.load(str(config_file))
#%%
para = gpu.Parameter(config, comm.get_number_of_process(), comm.get_pid())



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
para.set_use_AMD(use_AMD)

para.set_main_kernel("TurbulentViscosityCumulantK17CompChim" if para.get_use_AMD() else "CummulantK17CompChim")

para.set_SGS_constant(0.083)

def init_func(coord_x, coord_y, coord_z):
    return [
        0.0, 
        (u_star/kappa*np.log(max(coord_z/z_0,0)+1) + 2*np.sin(np.pi*16*coord_x/length[0])*np.sin(np.pi*8*coord_z/length[2]))/((coord_z/reference_height)**2+0.1)*dt/dx, 
        2*np.sin(np.pi*16*coord_x/length[0])*np.sin(np.pi*8*coord_z/length[2])/((coord_z/reference_height)**2+0.1)*dt/dx, 
        8*u_star/kappa*(np.sin(np.pi*8*coord_y/reference_height)*np.sin(np.pi*8*coord_z/reference_height)+np.sin(np.pi*8*coord_x/length[0]))/((length[2]/2-coord_z)**2+0.1)*dt/dx
        ]

para.set_initial_condition(init_func)
para.set_t_out(int(t_out/dt))
para.set_t_end(int(t_end/dt))
para.set_is_body_force(True)
para.set_has_wall_model_monitor(True)


grid_builder.add_coarse_grid(0.0, 0.0, 0.0, *length, dx)
grid_builder.set_periodic_boundary_condition(True, True, False)
grid_builder.build_grids(basics.LbmOrGks.LBM, False)
#%%
sampling_offset = 2
grid_builder.set_stress_boundary_condition(gpu.SideType.MZ, 0.0, 0.0, 1.0, sampling_offset, z_0/dx)
grid_builder.set_slip_boundary_condition(gpu.SideType.PZ, 0.0, 0.0, 0.0)

#%%
cuda_memory_manager = gpu.CudaMemoryManager(para)
grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager, comm)

#%%
wall_probe = gpu.probes.WallModelProbe("wallModelProbe", str(output_path), int(t_start_averaging/dt), int(t_start_tmp_averaging/dt), int(t_averaging/dt/4), int(t_start_out_probe/dt), int(t_out_probe/dt))
wall_probe.add_all_available_statistics()
wall_probe.set_file_name_to_n_out()
wall_probe.set_force_output_to_stress(True)
if para.get_is_body_force():
    wall_probe.set_evaluate_pressure_gradient(True)
planar_probe = gpu.probes.PlanarAverageProbe("planarAverageProbe", str(output_path), int(t_start_averaging/dt), int(t_start_tmp_averaging/dt), int(t_averaging/dt), int(t_start_out_probe/dt), int(t_out_probe/dt), "z")
para.add_probe(wall_probe)

#%%
sim = gpu.Simulation(para, cuda_memory_manager, comm, grid_generator)
#%%
sim.run()
MPI.Finalize()