from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBCAdapter, NoSlipBCAlgorithm
from pyfluids.geometry import GbCuboid3D, State
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import SimulationParameters, GridParameters, PhysicalParameters
from pyfluids.writer import Writer, WriterType

grid_params = GridParameters()
grid_params.delta_x = 1
grid_params.number_of_nodes_per_direction = [2, 2, 10]
grid_params.blocks_per_direction = [1, 1, 1]
grid_params.periodic_boundary_in_x1 = True
grid_params.periodic_boundary_in_x2 = True

physical_params = PhysicalParameters()
physical_params.lattice_viscosity = 0.005

sim_params = SimulationParameters()
sim_params.number_of_threads = 4
sim_params.number_of_timesteps = 10000
sim_params.timestep_log_interval = 1000


def run_simulation(physical_params=physical_params, grid_params=grid_params, sim_params=sim_params):
    simulation = Simulation()

    kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
    kernel.use_forcing = True
    kernel.forcing_in_x1 = 1e-6

    g_min_x1, g_min_x2, g_min_x3 = 0, 0, 0
    g_max_x1 = (grid_params.number_of_nodes_per_direction[0]) * grid_params.delta_x
    g_max_x2 = (grid_params.number_of_nodes_per_direction[1]) * grid_params.delta_x
    g_max_x3 = (grid_params.number_of_nodes_per_direction[2]) * grid_params.delta_x

    writer = Writer()
    writer.output_path = "./output"
    writer.type = WriterType.BINARY

    simulation.set_kernel_config(kernel)
    simulation.set_physical_parameters(physical_params)
    simulation.set_grid_parameters(grid_params)
    simulation.set_simulation_parameters(sim_params)
    simulation.set_writer(writer)

    no_slip_adapter = NoSlipBCAdapter()
    no_slip_adapter.algorithm = NoSlipBCAlgorithm()
    simulation.add_bc_adapter(no_slip_adapter)

    block_length = 3 * grid_params.delta_x
    simulation.add_object(
        GbCuboid3D(g_min_x1 - block_length,
                   g_min_x2 - block_length,
                   g_min_x3 - block_length,
                   g_max_x1 + block_length,
                   g_max_x2 + block_length,
                   g_min_x3),
        no_slip_adapter,
        State.SOLID, "/geo/addWallZMin")

    simulation.add_object(
        GbCuboid3D(g_min_x1 - block_length,
                   g_min_x2 - block_length,
                   g_max_x3,
                   g_max_x1 + block_length,
                   g_max_x2 + block_length,
                   g_max_x3 + block_length),
        no_slip_adapter,
        State.SOLID, "/geo/addWallZMax")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
