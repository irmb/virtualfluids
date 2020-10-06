from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBCAlgorithm, NoSlipBCAdapter, VelocityBCAdapter, DensityBCAdapter, \
    VelocityBCAlgorithm, NonReflectingOutflowBCAlgorithm
from pyfluids.geometry import GbCuboid3D
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import PhysicalParameters, SimulationParameters, GridParameters
from pyfluids.writer import Writer, WriterType
from pymuparser import Parser


def get_max_length(number_of_nodes_per_direction, delta_x):
    return (number_of_nodes_per_direction[0] * delta_x,
            number_of_nodes_per_direction[1] * delta_x,
            number_of_nodes_per_direction[2] * delta_x)


physical_params = PhysicalParameters()
physical_params.lattice_viscosity = 0.005

grid_params = GridParameters()
grid_params.number_of_nodes_per_direction = [200, 120, 120]
grid_params.blocks_per_direction = [2, 2, 2]
grid_params.delta_x = 0.125
grid_params.periodic_boundary_in_x1 = False
grid_params.periodic_boundary_in_x2 = True
grid_params.periodic_boundary_in_x3 = True

sim_params = SimulationParameters()
sim_params.timestep_log_interval = 1000
sim_params.number_of_timesteps = 1000
sim_params.number_of_threads = 4


def run_simulation(physical_parameters=physical_params, grid_parameters=grid_params,
                   sim_parameters=sim_params):
    wall_thickness = 3 * grid_parameters.delta_x

    min_x, min_y, min_z = 0, 0, 0
    max_x, max_y, max_z = get_max_length(grid_parameters.number_of_nodes_per_direction, grid_parameters.delta_x)

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

    no_slip_bc = NoSlipBCAdapter()
    no_slip_bc.algorithm = NoSlipBCAlgorithm()

    outflow_bc = DensityBCAdapter()
    outflow_bc.algorithm = NonReflectingOutflowBCAlgorithm()

    velocity_function = Parser()
    velocity_function.define_constant("u", 0.07)
    velocity_function.expression = "u"
    velocity_bc = VelocityBCAdapter(True, False, False, velocity_function, 0, -10)
    velocity_bc.algorithm = VelocityBCAlgorithm()

    kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)
    # kernel.use_forcing = True
    # kernel.forcing_in_x1 = 3e-6

    writer = Writer()
    writer.output_path = "./output"
    writer.type = WriterType.BINARY

    builder = Simulation()
    builder.set_writer(writer)

    builder.set_physical_parameters(physical_parameters)
    builder.set_grid_parameters(grid_parameters)
    builder.set_simulation_parameters(sim_parameters)
    builder.set_kernel_config(kernel)

    # builder.add_object(bottom_wall, no_slip_bc, 1, "/geo/bottomWall")
    # builder.add_object(top_wall, no_slip_bc, 1, "/geo/topWall")
    # builder.add_object(left_wall, no_slip_bc, 1, "/geo/leftWall")
    # builder.add_object(right_wall, no_slip_bc, 1, "/geo/rightWall")

    builder.add_object(obstacle, no_slip_bc, 1, "/geo/obstacle")

    builder.add_object(outflow_boundary, outflow_bc, 1, "/geo/outflow")
    builder.add_object(velocity_boundary, velocity_bc, 1, "/geo/velocityBoundary")

    builder.run_simulation()
