#include <pybind11/pybind11.h>
#include <gpu/GridGenerator/grid/BoundaryConditions/Side.h>

namespace boundary_conditions
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<SideType>(parentModule, "SideType")
        .value("MX", SideType::MX)
        .value("PX", SideType::PX)
        .value("MY", SideType::MY)
        .value("PY", SideType::PY)
        .value("MZ", SideType::MZ)
        .value("PZ", SideType::PZ)
        .value("GEOMETRY", SideType::GEOMETRY)
        .export_values();
    }
}