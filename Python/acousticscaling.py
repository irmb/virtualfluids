from pyfluids.kernel import LBMKernel
from pyfluids.parameters import GridParameters, PhysicalParameters, RuntimeParameters


class OneDirectionalAcousticScaling:

    def __init__(self, grid_parameters: GridParameters,
                 physical_parameters: PhysicalParameters,
                 runtime_parameters: RuntimeParameters,
                 kernel: LBMKernel):
        self._grid_params = grid_parameters
        self._physical_params = physical_parameters
        self._runtime_params = runtime_parameters
        self._kernel = kernel

    def configuration_for_scale_level(self, level: int = 1) -> (GridParameters,
                                                                PhysicalParameters,
                                                                RuntimeParameters,
                                                                LBMKernel):
        if level < 0:
            raise ValueError("level must be >= 0")

        grid_params = self.clone_grid_params_for_level(level)
        physical_params = self.clone_physical_parameters(level)
        runtime_params = self.clone_runtime_params_for_level(level)
        kernel = self.clone_kernel_for_level(level)

        return grid_params, physical_params, runtime_params, kernel

    def clone_grid_params_for_level(self, level) -> GridParameters:
        grid_params = GridParameters()
        grid_params.reference_direction_index = self._grid_params.reference_direction_index
        grid_params.periodic_boundary_in_x1 = self._grid_params.periodic_boundary_in_x1
        grid_params.periodic_boundary_in_x2 = self._grid_params.periodic_boundary_in_x2
        grid_params.periodic_boundary_in_x3 = self._grid_params.periodic_boundary_in_x3

        grid_params.number_of_nodes_per_direction = list(self._grid_params.number_of_nodes_per_direction)
        grid_params.blocks_per_direction = list(self._grid_params.blocks_per_direction)
        grid_params.node_distance = self._grid_params.node_distance

        if level > 0:
            grid_params.node_distance /= (level * 2)
            grid_params.number_of_nodes_per_direction = [grid_params.number_of_nodes_per_direction[0],
                                                         grid_params.number_of_nodes_per_direction[1],
                                                         grid_params.number_of_nodes_per_direction[2] * (level * 2)]

            grid_params.blocks_per_direction = [grid_params.blocks_per_direction[0],
                                                grid_params.blocks_per_direction[1],
                                                grid_params.blocks_per_direction[2] * (level * 2)]

        return grid_params

    def clone_physical_parameters(self, level):
        physical_params = PhysicalParameters()
        physical_params.lattice_viscosity = self._physical_params.lattice_viscosity

        if level > 0:
            physical_params.lattice_viscosity *= (level * 2)

        return physical_params

    def clone_runtime_params_for_level(self, level):
        runtime_params = RuntimeParameters()
        runtime_params.number_of_timesteps = self._runtime_params.number_of_timesteps
        runtime_params.number_of_threads = self._runtime_params.number_of_threads
        runtime_params.timestep_log_interval = self._runtime_params.timestep_log_interval

        if level > 0:
            runtime_params.number_of_timesteps *= (level * 2)

        return runtime_params

    def clone_kernel_for_level(self, level):
        kernel = LBMKernel(self._kernel.type)
        kernel.use_forcing = self._kernel.use_forcing
        kernel.forcing_in_x1 = self._kernel.forcing_in_x1
        kernel.forcing_in_x2 = self._kernel.forcing_in_x2
        kernel.forcing_in_x3 = self._kernel.forcing_in_x3

        if level > 0:
            kernel.forcing_in_x1 /= (level * 2)
            kernel.forcing_in_x2 /= (level * 2)
            kernel.forcing_in_x3 /= (level * 2)

        return kernel
