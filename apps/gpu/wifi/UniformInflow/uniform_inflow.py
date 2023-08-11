#%%
import numpy as np
from pathlib import Path
from mpi4py import MPI
from pyfluids import basics, gpu, logger
from wiFI.wind_farm import create_wind_farm_from_json
from wiFI.logging.logger import LoggerConfig
from wiFI.aeroelastics.stiff_rotor import StiffRotorGPU
from wiFI.interfaces.implementations.velocity_provider.VirtualFluids import VFFarm 
import multiprocessing as mp

mp.set_start_method("spawn", force=True)
#%%
def main():
    communicator = gpu.MpiCommunicator.get_instance()
    sim_name = "NREL5MW"
    sim_dir = Path("/workspaces/VirtualFluids_dev/output/wifi/")
    config_file = Path("apps/gpu/wifi/UniformInflow")/"configUniformInflow.txt"
    farm_file = Path("/workspaces/VirtualFluids_dev/wifi/resources/turbine_data/NREL5MW")/"SingleTurbine.json"
    use_tip_correction = False
    tip_speed_ratio = 7.5
    #%%
    logger.Logger.initialize_logger()
    #%%
    grid_builder = gpu.grid_generator.MultipleGridBuilder()

    config = basics.ConfigurationFile()
    config.load(str(config_file))

    para = gpu.Parameter(communicator.get_number_of_process(), communicator.get_pid(), config)
    para.set_use_streams(True)
    bc_factory = gpu.BoundaryConditionFactory()

    grid_scaling_factory = gpu.GridScalingFactory()
    grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)

    #%%
    turbine_diameter = config.get_float_value("turbineDiameter", 126)


    viscosity = config.get_float_value("viscosity", 1.56e-5)

    velocity  = 8
    mach = config.get_float_value("Ma", 0.05)
    nodes_per_diameter = config.get_uint_value("NodesPerDiameter", 32)

    density = config.get_float_value("Density", 1.225)
    level = 0
    n_blade_nodes  = config.get_int_value("NumberOfNodesPerAL", 32)


    # all in s
    t_start_out   = config.get_float_value("tStartOut")
    t_out         = config.get_float_value("tOut")
    t_end         = config.get_float_value("tEnd") # total time of simulation

    t_start_averaging      = config.get_float_value("tStartAveraging")
    t_start_tmp_averaging  = config.get_float_value("tStartTmpAveraging")
    t_averaging            = config.get_float_value("tAveraging")
    t_start_out_probe      = config.get_float_value("tStartOutProbe")
    t_out_probe            = config.get_float_value("tOutProbe")

    #%%
    length = np.array([4,3,3])*turbine_diameter
    dx = turbine_diameter / nodes_per_diameter
    dt = dx * mach / (np.sqrt(3) * velocity)
    velocity_LB = velocity * dt / dx # LB units
    viscosity_LB = viscosity * dt / (dx * dx) # LB units
    pressure_gradient = 0
    epsilon = dx*pow(2,-level)*1.5

    logger.vf_log_info(f"velocity  [dx/dt] = {velocity_LB}")
    logger.vf_log_info(f"dt   = {dt}")
    logger.vf_log_info(f"dx   = {dx}")
    logger.vf_log_info(f"viscosity [10^8 dx^2/dt] = {viscosity_LB*1e8}")
    logger.vf_log_info(f"dpdx  = {pressure_gradient}")
    logger.vf_log_info(f"mach number  = {mach}")

    farm = create_wind_farm_from_json(farm_file, sim_dir, tip_speed_ratio, velocity, True, log_turbine=True, logger_config=LoggerConfig(0, 1.0, timesteps_in_buffer=1))
    
    farm.turbine.add_blade_forces_logging(True)
    farm.turbine.add_blade_coordinate_logging(True)
    farm.turbine.add_blade_velocities_logging(True)

    #%%
    para.set_output_prefix(sim_name)
    para.set_print_files(True)

    para.set_forcing(0, 0, 0)
    para.set_velocity_LB(velocity_LB)
    para.set_viscosity_LB(viscosity_LB)    
    para.set_velocity_ratio(dx/dt)
    para.set_viscosity_ratio(dx*dx/dt)

    para.set_main_kernel("CumulantK17")

    para.set_timestep_start_out(int(t_start_out))
    # para.set_timestep_out(20)
    para.set_timestep_out(int(t_out))
    para.set_timestep_end(int(t_end))
    para.set_is_body_force(True)
    #%%
    tm_factory = gpu.TurbulenceModelFactory(para)
    tm_factory.read_config_file(config)
    #%%
    grid_builder.add_coarse_grid(-1.*turbine_diameter, -0.5 * length[1], -0.5 * length[2], length[0]-1.*turbine_diameter, 0.5 * length[1], 0.5 * length[2], dx)
    grid_builder.set_periodic_boundary_condition(False, True, True)
    grid_builder.build_grids(False)



    grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_LB, 0.0, 0.0)
    grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 0)

    bc_factory.set_velocity_boundary_condition(gpu.VelocityBC.VelocityCompressible)
    bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)

    #%%
    para.set_initial_condition_uniform(velocity_LB, 0, 0)

    coupled_farm = VFFarm(farm, density, epsilon, level, dt, dx,  n_blade_nodes, StiffRotorGPU, (density, ), use_tip_correction, 0)
    para.add_actuator(coupled_farm)

    # wall_model_probe = gpu.probes.WallModelProbe("wallModelProbe", para.get_output_path(), int(t_start_averaging/dt), int(t_start_tmp_averaging/dt), int(t_averaging/dt), int(t_start_out_probe/dt), int(t_out_probe/dt))
    # wall_model_probe.add_all_available_statistics()
    # wall_model_probe.set_file_name_to_n_out()
    # wall_model_probe.set_force_output_to_stress(True)
    # if para.get_is_body_force():
    #     wall_model_probe.set_evaluate_pressure_gradient(True)
    # para.add_probe(wall_model_probe)

    # plane_locs = [farm.positions.x[0] + i*turbine_diameter for i in range(-1,6)]

    # for n_probe, probe_pos in enumerate(plane_locs):
    #     plane_probe = gpu.probes.PlaneProbe(f"planeProbe_{n_probe+1}", para.get_output_path(), int(t_start_averaging/dt), int(t_averaging/dt), int(t_start_out_probe/dt), int(t_out_probe/dt))
    #     plane_probe.set_probe_plane(probe_pos, 0, 0, dx, length[1], length[2])
    #     plane_probe.add_all_available_statistics()
    #     para.add_probe(plane_probe)

    #%%
    cuda_memory_manager = gpu.CudaMemoryManager(para)
    grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager, communicator)
    #%%
    sim = gpu.Simulation(para, cuda_memory_manager, communicator, grid_generator, bc_factory, tm_factory, grid_scaling_factory)
    #%%
    sim.run()
    MPI.Finalize()

if __name__ == "__main__":
    main()