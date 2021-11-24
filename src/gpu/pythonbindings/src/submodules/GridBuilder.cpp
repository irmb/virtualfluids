#include <pybind11/pybind11.h>
#include <GridGenerator/grid/GridBuilder/MultipleGridBuilder.h>

namespace grid_builder{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {

        py::module gridBuilderModule = parentModule.def_submodule("grid_builder");

        py::class_<MultipleGridBuilder>(gridBuilderModule, "GridBuilder")
        .def(py::init<>())
        .def("add_coarse_grid", &MultipleGridBuilder::addCoarseGrid)
        .def("add_grid", &MultipleGridBuilder::addGrid)
        .def("add_geometry", &MultipleGridBuilder::addGeometry)
        .def("set_periodic_slip_conditions", &MultipleGridBuilder::setSlipBoundaryCondition)
        .def("set_periodic_velocity_conditions", &MultipleGridBuilder::setVelocityBoundaryCondition)
        .def("set_periodic_pressure_conditions", &MultipleGridBuilder::setPressureBoundaryCondition)
        .def("set_periodic_boundary_conditions", &MultipleGridBuilder::setPeriodicBoundaryCondition)
        .def("set_periodic_noslip_conditions", &MultipleGridBuilder::setNoSlipBoundaryCondition)
        .def("build_grids", &MultipleGridBuilder::buildGrids);
    }
}