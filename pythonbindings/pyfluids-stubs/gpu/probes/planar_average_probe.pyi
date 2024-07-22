from __future__ import annotations

import enum
from ... import gpu
from ...basics.geometry3d import Axis

class Statistic(enum.Enum):
    Means: int
    Covariances: int
    Skewness: int
    Flatness: int

class PlanarAverageProbe(gpu.Sampler):
    def __init__(
        self,
        para: gpu.Parameter,
        cuda_memory_manager: gpu.CudaMemoryManager,
        probe_name: str,
        output_path: str,
        t_start_avg: int,
        t_start_tmp_avg: int,
        t_avg: int,
        t_start_out: int,
        t_out: int,
        plane_normal: Axis,
        compute_time_averages: bool,
        compute_statistics_of_concentration: bool,
    ) -> None: ...
    def add_statistic(self, variable: Statistic) -> None: ...
    def add_all_available_statistics(self) -> None: ...
    def set_file_name_to_n_out(self) -> None: ...
