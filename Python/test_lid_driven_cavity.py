from pyfluids import Simulation
from pyfluids.boundaryconditions import NoSlipBCAdapter, NoSlipBCAlgorithm, VelocityBCAdapter, VelocityBCAlgorithm
from pyfluids.geometry import GbCuboid3D
from pyfluids.kernel import LBMKernel, KernelType
from pyfluids.parameters import GridParameters, PhysicalParameters, SimulationParameters
from pymuparser import Parser

builder = Simulation()
kernel = LBMKernel(KernelType.CompressibleCumulantFourthOrderViscosity)

sim_parameters = SimulationParameters()
sim_parameters.number_of_threads = 4
sim_parameters.number_of_timesteps = 10000
sim_parameters.timestep_log_interval = 1000

physical_parameters = PhysicalParameters()
physical_parameters.lattice_viscosity = 0.005

grid_parameters = GridParameters()
grid_parameters.number_of_nodes_per_direction = [64, 64, 64]
grid_parameters.blocks_per_direction = [2, 2, 2]
grid_parameters.delta_x = 1 / 10

builder.output_path = "./output"
builder.set_grid_parameters(grid_parameters)
builder.set_physical_parameters(physical_parameters)
builder.set_simulation_parameters(sim_parameters)
builder.set_kernel_config(kernel)

no_slip_bc_adapter = NoSlipBCAdapter()
no_slip_bc_adapter.algorithm = NoSlipBCAlgorithm()

fct = Parser()
fct.expression = "u"
fct.define_constant("u", 0.005)
velocity_bc_adapter = VelocityBCAdapter(True, True, False, fct, 0, -10.0)
velocity_bc_adapter.algorithm = VelocityBCAlgorithm()

g_minX1 = 0
g_minX2 = 0
g_minX3 = 0
g_maxX1 = grid_parameters.number_of_nodes_per_direction[0] * grid_parameters.delta_x
g_maxX2 = grid_parameters.number_of_nodes_per_direction[1] * grid_parameters.delta_x
g_maxX3 = grid_parameters.number_of_nodes_per_direction[2] * grid_parameters.delta_x

dx = grid_parameters.delta_x

wall_x_min = GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_minX1, g_maxX2 + dx, g_maxX3)
wall_x_max = GbCuboid3D(g_maxX1, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3)
wall_y_min = GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_minX2, g_maxX3)
wall_y_max = GbCuboid3D(g_minX1 - dx, g_maxX2, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_maxX3)
wall_z_min = GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_minX3 - dx, g_maxX1 + dx, g_maxX2 + dx, g_minX3)
wall_z_max = GbCuboid3D(g_minX1 - dx, g_minX2 - dx, g_maxX3, g_maxX1 + dx, g_maxX2 + dx, g_maxX3 + dx)

builder.add_object(wall_x_min, no_slip_bc_adapter, 1, "/geo/wallXmin")
builder.add_object(wall_x_max, no_slip_bc_adapter, 1, "/geo/wallXmax")
builder.add_object(wall_y_min, no_slip_bc_adapter, 1, "/geo/wallYmin")
builder.add_object(wall_y_max, no_slip_bc_adapter, 1, "/geo/wallYmax")
builder.add_object(wall_z_min, no_slip_bc_adapter, 1, "/geo/wallZmin")
builder.add_object(wall_z_max, velocity_bc_adapter, 1, "/geo/wallZmax")

builder.run_simulation()
