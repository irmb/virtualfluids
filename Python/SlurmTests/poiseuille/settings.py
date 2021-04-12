import os
from acousticscaling import AcousticScaling
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import RuntimeParameters, GridParameters, PhysicalParameters


grid_params = GridParameters()
grid_params.node_distance = 1
grid_params.number_of_nodes_per_direction = [1, 1, 16]
grid_params.blocks_per_direction = [1, 1, 8]
grid_params.periodic_boundary_in_x1 = True
grid_params.periodic_boundary_in_x2 = True

physical_params = PhysicalParameters()
physical_params.lattice_viscosity = 1e-3

runtime_params = RuntimeParameters()
runtime_params.number_of_threads = int(os.environ["SLURM_NTASKS"])
runtime_params.number_of_timesteps = 2_500_000
runtime_params.timestep_log_interval = 500_000

kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
kernel.use_forcing = True
kernel.forcing_in_x1 = 1e-9

Scaling = AcousticScaling(grid_params, physical_params, runtime_params, kernel)
