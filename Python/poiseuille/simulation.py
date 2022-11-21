from pyfluids import cpu


default_grid_params = cpu.parameters.GridParameters()
default_grid_params.node_distance = 10 / 32
default_grid_params.number_of_nodes_per_direction = [8, 8, 32]
default_grid_params.blocks_per_direction = [1, 1, 4]
default_grid_params.periodic_boundary_in_x1 = True
default_grid_params.periodic_boundary_in_x2 = True

default_physical_params = cpu.parameters.PhysicalParameters()
default_physical_params.lattice_viscosity = 0.005

default_runtime_params = cpu.parameters.RuntimeParameters()
default_runtime_params.number_of_threads = 4
default_runtime_params.number_of_timesteps = 10000
default_runtime_params.timestep_log_interval = 1000

default_kernel = cpu.kernel.LBMKernel(cpu.kernel.KernelType.CompressibleCumulantFourthOrderViscosity)
default_kernel.use_forcing = True
default_kernel.forcing_in_x1 = 1e-8

default_writer = cpu.writer.Writer()
default_writer.output_path = "./output"
default_writer.output_format = cpu.writer.OutputFormat.BINARY


default_kernel = cpu.kernel.LBMKernel(cpu.kernel.KernelType.CompressibleCumulantFourthOrderViscosity)
default_kernel.use_forcing = True
default_kernel.forcing_in_x1 = 1e-8


def run_simulation(physical_params=default_physical_params,
                   grid_params=default_grid_params,
                   runtime_params=default_runtime_params,
                   kernel=default_kernel,
                   writer=default_writer):
    simulation = cpu.Simulation()

    simulation.set_kernel_config(kernel)
    simulation.set_physical_parameters(physical_params)
    simulation.set_grid_parameters(grid_params)
    simulation.set_runtime_parameters(runtime_params)
    simulation.set_writer(writer)

    no_slip_bc = cpu.boundaryconditions.NoSlipBoundaryCondition()

    block_thickness = 3 * grid_params.node_distance
    simulation.add_object(
        cpu.geometry.GbCuboid3D(
            grid_params.bounding_box.min_x1 - block_thickness,
            grid_params.bounding_box.min_x2 - block_thickness,
            grid_params.bounding_box.min_x3 - block_thickness,
            grid_params.bounding_box.max_x1 + block_thickness,
            grid_params.bounding_box.max_x2 + block_thickness,
            grid_params.bounding_box.min_x3),
        no_slip_bc,
        cpu.geometry.State.SOLID, "/geo/addWallZMin")

    simulation.add_object(
        cpu.geometry.GbCuboid3D(
            grid_params.bounding_box.min_x1 - block_thickness,
            grid_params.bounding_box.min_x2 - block_thickness,
            grid_params.bounding_box.max_x3,
            grid_params.bounding_box.max_x1 + block_thickness,
            grid_params.bounding_box.max_x2 + block_thickness,
            grid_params.bounding_box.max_x3 + block_thickness),
        no_slip_bc,
        cpu.geometry.State.SOLID, "/geo/addWallZMax")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
