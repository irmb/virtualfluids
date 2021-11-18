import numpy as np
import scipy.stats as stats
import errors
from SlurmTests.poiseuille.result_collector import collect_results
from SlurmTests.poiseuille.settings import Scaling

analytical_results, numerical_results = collect_results()
normalized_l2_errors = [errors.normalized_l2_error(analytical, numerical)
                        for analytical, numerical in zip(analytical_results, numerical_results)]

nodes_in_x3_per_run = []
for simulation_run in range(0, 3):
    grid_params, _, _, _ = Scaling.configuration_for_scale_level(simulation_run)
    nodes_in_x3_per_run.append(grid_params.number_of_nodes_per_direction[2])

nodes_as_log = [np.log10(node) for node in nodes_in_x3_per_run]
l2_norms_as_log = [np.log10(l2) for l2 in normalized_l2_errors]
res = stats.linregress(nodes_as_log, l2_norms_as_log)

assert res.slope <= -2, f"Expected slope of l2 error to be <= -2, but was {res.slope}"
