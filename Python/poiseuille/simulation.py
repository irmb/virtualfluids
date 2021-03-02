from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBoundaryCondition
from pyfluids.geometry import GbCuboid3D, State
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import RuntimeParameters, GridParameters, PhysicalParameters
from pyfluids.writer import Writer, OutputFormat

default_grid_params = GridParameters()
default_grid_params.node_distance = 10 / 32
default_grid_params.number_of_nodes_per_direction = [8, 8, 32]
default_grid_params.blocks_per_direction = [1, 1, 4]
default_grid_params.periodic_boundary_in_x1 = True
default_grid_params.periodic_boundary_in_x2 = True

default_physical_params = PhysicalParameters()
default_physical_params.lattice_viscosity = 0.005

default_runtime_params = RuntimeParameters()
default_runtime_params.number_of_threads = 4
default_runtime_params.number_of_timesteps = 10000
default_runtime_params.timestep_log_interval = 1000


default_kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
default_kernel.use_forcing = True
default_kernel.forcing_in_x1 = 1e-8


def run_simulation(physical_params=default_physical_params,
                   grid_params=default_grid_params,
                   runtime_params=default_runtime_params,
                   kernel=default_kernel):
    simulation = Simulation()

    writer = Writer()
    writer.output_path = "./output"
    writer.output_format = OutputFormat.BINARY

    simulation.set_kernel_config(kernel)
    simulation.set_physical_parameters(physical_params)
    simulation.set_grid_parameters(grid_params)
    simulation.set_runtime_parameters(runtime_params)
    simulation.set_writer(writer)

    no_slip_bc = NoSlipBoundaryCondition()

    block_thickness = 3 * grid_params.node_distance
    simulation.add_object(
        GbCuboid3D(
            grid_params.bounding_box.min_x1 - block_thickness,
            grid_params.bounding_box.min_x2 - block_thickness,
            grid_params.bounding_box.min_x3 - block_thickness,
            grid_params.bounding_box.max_x1 + block_thickness,
            grid_params.bounding_box.max_x2 + block_thickness,
            grid_params.bounding_box.min_x3),
        no_slip_bc,
        State.SOLID, "/geo/addWallZMin")

    simulation.add_object(
        GbCuboid3D(
            grid_params.bounding_box.min_x1 - block_thickness,
            grid_params.bounding_box.min_x2 - block_thickness,
            grid_params.bounding_box.max_x3,
            grid_params.bounding_box.max_x1 + block_thickness,
            grid_params.bounding_box.max_x2 + block_thickness,
            grid_params.bounding_box.max_x3 + block_thickness),
        no_slip_bc,
        State.SOLID, "/geo/addWallZMax")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
