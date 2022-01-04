#include <pybind11/pybind11.h>
#include "submodules/actuator_line.cpp"
#include "submodules/simulation.cpp"
#include "submodules/parameter.cpp"
#include "submodules/boundary_conditions.cpp"
#include "submodules/communicator.cpp"
#include "submodules/grid_builder.cpp"

namespace gpu
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module gpuModule = parentModule.def_submodule("gpu");
        simulation::makeModule(gpuModule);
        parameter::makeModule(gpuModule);
        actuator_line::makeModule(gpuModule);
        boundary_conditions::makeModule(gpuModule);
        communicator::makeModule(gpuModule);
        grid_builder::makeModule(gpuModule);

        return gpuModule;
    }
}