from pyfluids.cpu import Simulation
from pyfluids.cpu.boundaryconditions import NoSlipBoundaryCondition, VelocityBoundaryCondition, DensityBoundaryCondition
from pyfluids.cpu.geometry import GbCuboid3D
from pyfluids.cpu.kernel import LBMKernel, KernelType
from pyfluids.cpu.parameters import PhysicalParameters, RuntimeParameters, GridParameters
from pyfluids.cpu.writer import Writer, OutputFormat
from pymuparser import Parser

import os


def get_max_length(number_of_nodes_per_direction, delta_x):
    return (number_of_nodes_per_direction[0] * delta_x,
            number_of_nodes_per_direction[1] * delta_x,
            number_of_nodes_per_direction[2] * delta_x)


physical_params = PhysicalParameters()
physical_params.lattice_viscosity = 0.005

grid_params = GridParameters()
grid_params.number_of_nodes_per_direction = [200, 120, 120]
grid_params.blocks_per_direction = [2, 2, 2]
grid_params.node_distance = 0.125
grid_params.periodic_boundary_in_x1 = False
grid_params.periodic_boundary_in_x2 = True
grid_params.periodic_boundary_in_x3 = True

runtime_params = RuntimeParameters()
runtime_params.timestep_log_interval = 1000
runtime_params.number_of_timesteps = 50000
runtime_params.number_of_threads = int(os.environ.get("OMP_NUM_THREADS", 4))


def run_simulation(physical_parameters=physical_params, grid_parameters=grid_params,
                   runtime_parameters=runtime_params):
    wall_thickness = 3 * grid_parameters.node_distance

    min_x, min_y, min_z = 0, 0, 0
    max_x, max_y, max_z = get_max_length(grid_parameters.number_of_nodes_per_direction, grid_parameters.node_distance)

    bottom_wall = GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, min_z, max_x + wall_thickness,
                             max_y + wall_thickness, min_z - wall_thickness)

    top_wall = GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, max_z, max_x + wall_thickness,
                          max_y + wall_thickness,
                          max_z + wall_thickness)

    left_wall = GbCuboid3D(min_x - wall_thickness, min_y, min_z - wall_thickness, max_x + wall_thickness,
                           min_y - wall_thickness,
                           max_z + wall_thickness)

    right_wall = GbCuboid3D(min_x - wall_thickness, max_y, min_z - wall_thickness, max_x + wall_thickness,
                            max_y + wall_thickness, max_z + wall_thickness)

    obstacle = GbCuboid3D(7, 7, 7, 8, 8, 8)

    velocity_boundary = GbCuboid3D(min_x - wall_thickness, min_y - wall_thickness, min_z - wall_thickness, min_x,
                                   max_y + wall_thickness, max_z + wall_thickness)

    outflow_boundary = GbCuboid3D(max_x, min_y - wall_thickness, min_z - wall_thickness, max_x + wall_thickness,
                                  max_y + wall_thickness, max_z + wall_thickness)

    no_slip_bc = NoSlipBoundaryCondition()

    outflow_bc = DensityBoundaryCondition()

    velocity_function = Parser()
    velocity_function.define_constant("u", 0.07)
    velocity_function.expression = "u"
    velocity_bc = VelocityBoundaryCondition(True, False, False, velocity_function, 0, -10)

    kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
    # kernel.use_forcing = True
    # kernel.forcing_in_x1 = 3e-6

    writer = Writer()
    writer.output_path = "./output"
    writer.output_format = OutputFormat.BINARY

    simulation = Simulation()
    simulation.set_writer(writer)

    simulation.set_physical_parameters(physical_parameters)
    simulation.set_grid_parameters(grid_parameters)
    simulation.set_runtime_parameters(runtime_parameters)
    simulation.set_kernel_config(kernel)

    # simulation.add_object(bottom_wall, no_slip_bc, 1, "/geo/bottomWall")
    # simulation.add_object(top_wall, no_slip_bc, 1, "/geo/topWall")
    # simulation.add_object(left_wall, no_slip_bc, 1, "/geo/leftWall")
    # simulation.add_object(right_wall, no_slip_bc, 1, "/geo/rightWall")

    simulation.add_object(obstacle, no_slip_bc, 1, "/geo/obstacle")

    simulation.add_object(outflow_boundary, outflow_bc, 1, "/geo/outflow")
    simulation.add_object(velocity_boundary, velocity_bc, 1, "/geo/velocityBoundary")

    simulation.run_simulation()


if __name__ == "__main__":
    run_simulation()
