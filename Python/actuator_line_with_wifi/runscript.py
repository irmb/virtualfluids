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
from pathlib import Path
from wiFI.turbine import create_turbine_from_farm_json
from wiFI.logging import LoggerConfig, LogEntry
from wiFI.interfaces.implementations.velocity_provider.VirtualFluids.VirtualFluids import create_standard_actuator_farm
from wiFI.controller.controller import ControllerTypes
from pyfluids import basics, gpu, logger, parallel
import numpy as np


#%%
sim_name = "ABL"
config_file = Path(__file__).parent/"configActuatorLine.txt"
turbine_file = Path(__file__).parent/"SingleTurbine.json"
controller_file = Path(__file__).parent/"controller.json"

def load_config_file(config_file: Path):
    config = basics.ConfigurationFile()
    config.load(str(config_file))
    return config

def add_plane_probes(plane_locs: np.ndarray, para: gpu.Parameter, t_start_averaging: float, t_start_out_probe: float, t_out_probe: float, dt: float, dx: float, length: np.ndarray):
    for n_probe, probe_pos in enumerate(plane_locs):
        plane_probe = gpu.probes.PlaneProbe(f"planeProbe_{n_probe+1}", para.get_output_path(), int(t_start_averaging/dt), 10, int(t_start_out_probe/dt), int(t_out_probe/dt))
        plane_probe.set_probe_plane(probe_pos, -length[1]/2, -length[2]/2, dx, length[1], length[2])
        plane_probe.add_all_available_statistics()
        para.add_probe(plane_probe)


def main(sim_name: str, config_file: Path, turbine_file: Path, controller_file: Path):
    logger.Logger.initialize_logger()

    config = load_config_file(config_file)
    grid_builder = gpu.grid_generator.MultipleGridBuilder()
    communicator = parallel.MPICommunicator.get_instance()

    para = gpu.Parameter(communicator.get_number_of_processes(), communicator.get_process_id(), config)
    bc_factory = gpu.BoundaryConditionFactory()


    viscosity = config.get_float_value("viscosity", 1.56e-5)
    velocity  = config.get_float_value("velocity", 9.0)

    mach = config.get_float_value("Ma", 0.1)
    nodes_per_diameter = config.get_uint_value("nodesPerDiameter", 32)


    density = config.get_float_value("Density", 1.225)
    level = 0
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


    turbine_model = create_turbine_from_farm_json(turbine_file,
                                                  n_blade_nodes,
                                                  load_data = True,
                                                  use_gpu = True)
    
    diameter =  turbine_model.blade_tip_radius*2

    length = np.array([24,8,8])*diameter
    dx = diameter/nodes_per_diameter
    dt = dx * mach / (np.sqrt(3) * velocity)
    velocity_ratio = dx/dt
    velocity_LB = velocity / velocity_ratio  # LB units
    viscosity_LB = viscosity / (velocity_ratio * dx)  # LB units

    logger.vf_log_info(f"velocity  [dx/dt] = {velocity_LB}")
    logger.vf_log_info(f"dt   = {dt}")
    logger.vf_log_info(f"dx   = {dx}")
    logger.vf_log_info(f"viscosity [10^8 dx^2/dt] = {viscosity_LB*1e8}")

    para.set_output_prefix(sim_name)
    para.set_print_files(True)
    output_path = Path(para.get_output_path())
    output_path.mkdir(exist_ok=True)

    para.set_forcing(0, 0, 0)
    para.set_velocity_LB(velocity_LB)
    para.set_viscosity_LB(viscosity_LB)    
    para.set_velocity_ratio(dx/dt)
    para.set_viscosity_ratio(dx*dx/dt)
    para.set_density_ratio(1.0)

    para.configure_main_kernel(gpu.kernel.compressible.K17CompressibleNavierStokes)

    para.set_timestep_start_out(int(t_start_out/dt))
    para.set_timestep_out(int(t_out/dt))
    para.set_timestep_end(int(t_end/dt))
    para.set_is_body_force(True)

    tm_factory = gpu.TurbulenceModelFactory(para)
    tm_factory.read_config_file(config)
    #%%
    grid_scaling_factory = gpu.GridScalingFactory()
    grid_scaling_factory.set_scaling_factory(gpu.GridScaling.ScaleCompressible)

    grid_builder.add_coarse_grid(-3*diameter, -length[1]/2, -length[2]/2, length[0]-3*diameter, length[1]/2, length[2]/2, dx)
    grid_builder.set_periodic_boundary_condition(False, True, True)
    grid_builder.build_grids(False)

    if read_precursor:
        precursor = gpu.create_file_collection(precursor_directory + "/precursor", gpu.FileType.VTK)
        grid_builder.set_precursor_boundary_condition(gpu.SideType.MX, precursor, nTReadPrecursor, 0, 0, 0)
    else:
        grid_builder.set_velocity_boundary_condition(gpu.SideType.MX, velocity_LB, 0, 0)

    grid_builder.set_pressure_boundary_condition(gpu.SideType.PX, 0)

    bc_factory.set_stress_boundary_condition(gpu.StressBC.StressPressureBounceBack)
    bc_factory.set_slip_boundary_condition(gpu.SlipBC.SlipBounceBack) 
    bc_factory.set_pressure_boundary_condition(gpu.PressureBC.OutflowNonReflective)
    if read_precursor:
        bc_factory.set_precursor_boundary_condition(gpu.PrecursorBC.DistributionsPrecursor if use_distributions else gpu.PrecursorBC.VelocityPrecursor)
    else:
        bc_factory.set_velocity_boundary_condition(gpu.VelocityBC.VelocityAndPressureCompressible)

    para.set_outflow_pressure_correction_factor(0.0)

    para.set_initial_condition_uniform(velocity_LB, 0, 0)


    logging_config = LoggerConfig("wifi", output_path, start_time=1000, log_period=100)
    logging_dict = {"wind_farm": [LogEntry("rotor_speed", True, True),
                                  LogEntry("azimuth", True, True)]}


    tip_speed_ratio = 7.55*np.ones(turbine_model.n_turbines)
    rotor_speeds = tip_speed_ratio * velocity / turbine_model.blade_tip_radius
    smearing_width = 2*dx
    farm = create_standard_actuator_farm(logging_config, logging_dict, turbine_model, density, smearing_width,
                                         level, dt, dx, communicator.get_process_id(),
                                         ControllerTypes.Greedy, controller_file, rotor_speeds)

    # farm = gpu.ActuatorFarmStandalone(turbine_model.blade_tip_radius*2, turbine_model.n_nodes_per_blade, turbine_model.hub_positions.x, turbine_model.hub_positions.y, turbine_model.hub_positions.z, rotor_speeds, density, smearing_width, level, dt, dx)
    farm.enable_output("ALM", 0, int(t_out_probe/dt))
    para.add_actuator(farm)
    
    plane_locs = np.array([-1,1,2,3,4])*diameter
    add_plane_probes(plane_locs, para, t_start_averaging, t_start_out_probe, t_out_probe, dt, dx, length)
    plane_probe = gpu.probes.PlaneProbe(f"streamwiseProbe", para.get_output_path(), int(t_start_averaging/dt), 10, int(t_start_out_probe/dt), int(t_out_probe/dt))
    plane_probe.set_probe_plane(-diameter*3, 0, -length[2]/2, length[0], dx, length[2])
    plane_probe.add_all_available_statistics()
    para.add_probe(plane_probe)

    cuda_memory_manager = gpu.CudaMemoryManager(para)
    grid_generator = gpu.GridProvider.make_grid_generator(grid_builder, para, cuda_memory_manager, communicator)
    sim = gpu.Simulation(para, cuda_memory_manager, communicator, grid_generator, bc_factory, tm_factory, grid_scaling_factory)
    sim.run()


if __name__ == '__main__':
    main(sim_name, config_file, turbine_file, controller_file)
