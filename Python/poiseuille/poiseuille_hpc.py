from poiseuille.simulation import run_simulation
from pyfluids.parameters import *

grid_parameters = GridParameters()
grid_parameters.number_of_nodes_per_direction = [32, 32, 128]
grid_parameters.node_distance = 1
grid_parameters.blocks_per_direction = [1, 1, 2]

physical_parameters = PhysicalParameters()
physical_parameters.lattice_viscosity = 0.0005

runtime_parameters = RuntimeParameters()
runtime_parameters.number_of_threads = 1
runtime_parameters.number_of_timesteps = 1000
runtime_parameters.timestep_log_interval = 100

run_simulation(physical_parameters, grid_parameters, runtime_parameters)