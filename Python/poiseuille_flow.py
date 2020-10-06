from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBCAdapter, NoSlipBCAlgorithm
from pyfluids.geometry import GbCuboid3D
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.writer import Writer, WriterType


def simulate_poiseuille_flow(physical_params, grid_params, sim_params):
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
        1, "/geo/addWallZMin")

    simulation.add_object(
        GbCuboid3D(g_min_x1 - block_length,
                   g_min_x2 - block_length,
                   g_max_x3,
                   g_max_x1 + block_length,
                   g_max_x2 + block_length,
                   g_max_x3 + block_length),
        no_slip_adapter,
        1, "/geo/addWallZMax")

    simulation.run_simulation()
