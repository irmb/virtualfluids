from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBoundaryCondition
from pyfluids.geometry import GbCuboid3D, State
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import RuntimeParameters, GridParameters, PhysicalParameters
from pyfluids.writer import Writer, OutputFormat

grid_params = GridParameters()
grid_params.node_distance = 1
grid_params.number_of_nodes_per_direction = [1, 1, 10]
grid_params.blocks_per_direction = [1, 1, 1]
grid_params.periodic_boundary_in_x1 = True
grid_params.periodic_boundary_in_x2 = True

physical_params = PhysicalParameters()
physical_params.lattice_viscosity = 0.005

runtime_params = RuntimeParameters()
runtime_params.number_of_threads = 4
runtime_params.number_of_timesteps = 10000
runtime_params.timestep_log_interval = 1000


def run_simulation(physical_params=physical_params, grid_params=grid_params, runtime_params=runtime_params):
    simulation = Simulation()

    kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
    kernel.use_forcing = True
    kernel.forcing_in_x1 = 1e-6

    node_distance = grid_params.node_distance
    min_x1, min_x2, min_x3 = 0, 0, 0
    max_x1 = grid_params.number_of_nodes_per_direction[0] * node_distance
    max_x2 = grid_params.number_of_nodes_per_direction[1] * node_distance
    max_x3 = grid_params.number_of_nodes_per_direction[2] * node_distance

    writer = Writer()
    writer.output_path = "./output"
    writer.output_format = OutputFormat.BINARY

    simulation.set_kernel_config(kernel)
    simulation.set_physical_parameters(physical_params)
    simulation.set_grid_parameters(grid_params)
    simulation.set_runtime_parameters(runtime_params)
    simulation.set_writer(writer)

    no_slip_bc = NoSlipBoundaryCondition()

    block_width = 3 * node_distance
    simulation.add_object(
        GbCuboid3D(min_x1 - block_width,
                   min_x2 - block_width,
                   min_x3 - block_width,
                   max_x1 + block_width,
                   max_x2 + block_width,
                   min_x3),
        no_slip_bc,
        State.SOLID, "/geo/addWallZMin")

    simulation.add_object(
        GbCuboid3D(min_x1 - block_width,
                   min_x2 - block_width,
                   max_x3,
                   max_x1 + block_width,
                   max_x2 + block_width,
                   max_x3 + block_width),
        no_slip_bc,
        State.SOLID, "/geo/addWallZMax")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
