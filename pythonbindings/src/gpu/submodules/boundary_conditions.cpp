#include <pybind11/pybind11.h>
#include <gpu/GridGenerator/grid/BoundaryConditions/Side.h>

namespace boundary_conditions
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module boundaryConditionsModule = parentModule.def_submodule("boundary_conditions");

        py::enum_<SideType>(boundaryConditionsModule, "side_type")
        .value("MX", SideType::MX)
        .value("PX", SideType::PX)
        .value("MY", SideType::MY)
        .value("PY", SideType::PY)
        .value("MZ", SideType::MZ)
        .value("PZ", SideType::PZ)
        .value("GEOMETRY", SideType::GEOMETRY)
        .export_values();

        return boundaryConditionsModule;
    }
}