from poiseuille.simulation import run_simulation
from pyfluids import cpu

grid_parameters = cpu.prameters.GridParameters()
grid_parameters.number_of_nodes_per_direction = [64, 64, 512]
grid_parameters.node_distance = 1
grid_parameters.blocks_per_direction = [1, 2, 2]

physical_parameters = cpu.prameters.PhysicalParameters()
physical_parameters.lattice_viscosity = 0.0005

runtime_parameters = cpu.prameters.RuntimeParameters()
runtime_parameters.number_of_threads = 4
runtime_parameters.number_of_timesteps = 1000
runtime_parameters.timestep_log_interval = 100

run_simulation(physical_parameters, grid_parameters, runtime_parameters)