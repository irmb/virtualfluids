import os

from SlurmTests.poiseuille.settings import Scaling
from poiseuille.simulation import run_simulation
from pyfluids.writer import Writer, OutputFormat


scale_level = int(os.environ["PYFLUIDS_SCALE_LEVEL"])
grid_params, physical_params, runtime_params, kernel = Scaling.configuration_for_scale_level(scale_level)

writer = Writer()
writer.output_format = OutputFormat.BINARY
writer.output_path = "./output-" + str(scale_level)

run_simulation(grid_params=grid_params,
               physical_params=physical_params,
               runtime_params=runtime_params,
               kernel=kernel,
               writer=writer)
