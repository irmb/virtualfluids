#include <pybind11/pybind11.h>
#include "submodules/actuator_line.cpp"
#include "submodules/pre_collision_interactor.cpp"
#include "submodules/simulation.cpp"
#include "submodules/parameter.cpp"
#include "submodules/boundary_conditions.cpp"
#include "submodules/communicator.cpp"
#include "submodules/grid_builder.cpp"
#include "submodules/cuda_memory_manager.cpp"
#include "submodules/grid_provider.cpp"
#include "submodules/probes.cpp"
#include "submodules/kernel_factory.cpp"
#include "submodules/pre_processor_factory.cpp"
#include "submodules/file_writer.cpp"

namespace gpu
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module gpuModule = parentModule.def_submodule("gpu");
        simulation::makeModule(gpuModule);
        parameter::makeModule(gpuModule);
        pre_collision_interactor::makeModule(gpuModule);
        actuator_line::makeModule(gpuModule);
        boundary_conditions::makeModule(gpuModule);
        communicator::makeModule(gpuModule); 
        grid_builder::makeModule(gpuModule);
        cuda_memory_manager::makeModule(gpuModule);
        grid_provider::makeModule(gpuModule);
        probes::makeModule(gpuModule);
        kernel_factory::makeModule(gpuModule);
        pre_processor_factory::makeModule(gpuModule);
        file_writer::makeModule(gpuModule);
        return gpuModule;
    }
}