#include <pybind11/pybind11.h>
#include "submodules/boundaryconditions.cpp"
#include "submodules/simulationconfig.cpp"
#include "submodules/geometry.cpp"
#include "submodules/kernel.cpp"
#include "submodules/simulationparameters.cpp"
#include "submodules/writer.cpp"

namespace cpu
{
    namespace py = pybind11;
    py::module makeModule(py::module_ &parentModule)
    {
        py::module cpuModule = parentModule.def_submodule("cpu");
        boundaryconditions::makeModule(cpuModule);
        simulation::makeModule(cpuModule);
        geometry::makeModule(cpuModule);
        kernel::makeModule(cpuModule);
        parameters::makeModule(cpuModule);
        writer::makeModule(cpuModule);
        return cpuModule;
    }
}