#include <pybind11/pybind11.h>
#include "submodules/actuator_line.cpp"
#include "submodules/pre_collision_interactor.cpp"
#include "submodules/simulation.cpp"
#include "submodules/parameter.cpp"
#include "submodules/boundary_conditions.cpp"
#include "submodules/communicator.cpp"
#include "submodules/cuda_memory_manager.cpp"
#include "submodules/probes.cpp"
#include "submodules/precursor_writer.cpp"
#include "submodules/grid_provider.cpp"
#include "submodules/grid_generator.cpp"
#include "submodules/turbulence_models.cpp"
#include "submodules/velocity_setter.cpp"

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
        velocity_setter::makeModule(gpuModule);
        communicator::makeModule(gpuModule); 
        cuda_memory_manager::makeModule(gpuModule);
        probes::makeModule(gpuModule);
        precursor_writer::makeModule(gpuModule);
        grid_generator::makeModule(gpuModule);
        grid_provider::makeModule(gpuModule);
        turbulence_model::makeModule(gpuModule);
        return gpuModule;
    }
}