#include <pybind11/pybind11.h>
#include "gpu/GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"
#include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/geometries/Object.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/Side.h"

namespace grid_builder
{

    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {        
        py::class_<GridBuilder, std::shared_ptr<GridBuilder>>(parentModule, "GridBuilder")
        .def("get_number_of_grid_levels", &GridBuilder::getNumberOfGridLevels)
        .def("get_grid", &GridBuilder::getGrid);

        py::class_<LevelGridBuilder, GridBuilder, std::shared_ptr<LevelGridBuilder>>(parentModule, "LevelGridBuilder")
        .def("get_grid", py::overload_cast<int, int>(&LevelGridBuilder::getGrid))
        .def("set_slip_boundary_condition", &LevelGridBuilder::setSlipBoundaryCondition)
        .def("set_velocity_boundary_condition", &LevelGridBuilder::setVelocityBoundaryCondition)
        .def("set_pressure_boundary_condition", &LevelGridBuilder::setPressureBoundaryCondition)
        .def("set_periodic_boundary_condition", &LevelGridBuilder::setPeriodicBoundaryCondition)
        .def("set_no_slip_boundary_condition", &LevelGridBuilder::setNoSlipBoundaryCondition);

        py::class_<MultipleGridBuilder, LevelGridBuilder, std::shared_ptr<MultipleGridBuilder>>(parentModule, "MultipleGridBuilder")
        .def("make_shared", &MultipleGridBuilder::makeShared, py::return_value_policy::reference)
        .def("add_coarse_grid", &MultipleGridBuilder::addCoarseGrid)
        .def("add_grid", py::overload_cast<Object*>(&MultipleGridBuilder::addGrid))
        .def("add_grid", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGrid))
        .def("add_geometry", py::overload_cast<Object*>(&MultipleGridBuilder::addGeometry))
        .def("add_geometry", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGeometry))
        .def("get_number_of_levels", &MultipleGridBuilder::getNumberOfLevels)
        .def("build_grids", &MultipleGridBuilder::buildGrids);
    }
}