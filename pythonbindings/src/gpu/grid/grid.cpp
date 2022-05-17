#include <pybind11/pybind11.h>
#include "grid_builder.cpp"
#include "grid_provider.cpp"
#include "grid_factory.cpp"
#include "grid_generator.cpp"

namespace grid
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module gridModule = parentModule.def_submodule("grid");

        grid_builder::makeModule(gridModule);
        grid_provider::makeModule(gridModule);
        grid_generator::makeModule(gridModule);
        grid_factory::makeModule(gridModule);
        return gridModule;
    }
}