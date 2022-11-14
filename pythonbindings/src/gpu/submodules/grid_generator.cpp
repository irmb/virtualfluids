#include <pybind11/pybind11.h>
#include "gpu/GridGenerator/utilities/communication.h"
#include "gpu/GridGenerator/geometries/Object.h"
#include "gpu/GridGenerator/geometries/BoundingBox/BoundingBox.h"
#include "gpu/GridGenerator/geometries/Conglomerate/Conglomerate.h"
#include "gpu/GridGenerator/geometries/Cuboid/Cuboid.h"
#include "gpu/GridGenerator/geometries/Sphere/Sphere.h"
#include "gpu/GridGenerator/geometries/TriangularMesh/TriangularMesh.h"
#include "gpu/GridGenerator/grid/GridFactory.h"
#include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"
#include "gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h"
#include "gpu/GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

namespace grid_generator
{
    namespace py = pybind11;
    py::module makeModule(py::module_ &parentModule)
    {  
        py::module gridGeneratorModule = parentModule.def_submodule("grid_generator");

        //TODO:
        // py::enum_<CommunicationDirections>(gridGeneratorModule, "CommunicationDirections")
        // .value("MX", CommunicationDirections::MX)
        // .value("PX", CommunicationDirections::PX)
        // .value("MY", CommunicationDirections::MY)
        // .value("PY", CommunicationDirections::PY)
        // .value("MZ", CommunicationDirections::MZ)
        // .value("PZ", CommunicationDirections::PZ);

        py::class_<GridFactory, std::shared_ptr<GridFactory>>(gridGeneratorModule, "GridFactory")
        .def("make", &GridFactory::make, py::return_value_policy::reference);

        py::class_<BoundingBox, std::shared_ptr<BoundingBox>>(gridGeneratorModule, "BoundingBox")
        .def(py::init<real, real, real, real, real, real>(),"min_x","max_x","min_y","max_y","min_z","max_z");

        py::class_<Object, std::shared_ptr<Object>>(gridGeneratorModule, "Object");
        
        py::class_<Conglomerate, Object, std::shared_ptr<Conglomerate>>(gridGeneratorModule, "Conglomerate")
        .def("make_shared", &Conglomerate::makeShared, py::return_value_policy::reference)
        .def("add", &Conglomerate::add)
        .def("subtract", &Conglomerate::subtract);

        py::class_<Cuboid, Object, std::shared_ptr<Cuboid>>(gridGeneratorModule, "Cuboid")
        .def(py::init<const double&, const double&, const double&, const double&, const double&, const double&>(),
                        "min_x1", "min_x2", "min_x3", "max_x1", "max_x2", "max_x3");

        py::class_<Sphere, Object, std::shared_ptr<Sphere>>(gridGeneratorModule, "Sphere")
        .def("make_shared", &Sphere::makeShared, py::return_value_policy::reference);

        py::class_<TriangularMesh, Object, std::shared_ptr<TriangularMesh>>(gridGeneratorModule, "TriangularMesh")
        .def("make", &TriangularMesh::make, py::return_value_policy::reference);

        py::class_<GridBuilder, std::shared_ptr<GridBuilder>>(gridGeneratorModule, "GridBuilder")
        .def("get_number_of_grid_levels", &GridBuilder::getNumberOfGridLevels)
        .def("get_grid", &GridBuilder::getGrid);

        py::class_<LevelGridBuilder, GridBuilder, std::shared_ptr<LevelGridBuilder>>(gridGeneratorModule, "LevelGridBuilder")
        .def("get_grid", py::overload_cast<int, int>(&LevelGridBuilder::getGrid))
        .def("set_slip_boundary_condition", &LevelGridBuilder::setSlipBoundaryCondition)
        .def("set_velocity_boundary_condition", &LevelGridBuilder::setVelocityBoundaryCondition)
        .def("set_pressure_boundary_condition", &LevelGridBuilder::setPressureBoundaryCondition)
        .def("set_periodic_boundary_condition", &LevelGridBuilder::setPeriodicBoundaryCondition)
        .def("set_no_slip_boundary_condition", &LevelGridBuilder::setNoSlipBoundaryCondition)
        .def("set_precursor_boundary_condition", &LevelGridBuilder::setPrecursorBoundaryCondition)
        .def("set_stress_boundary_condition", &LevelGridBuilder::setStressBoundaryCondition);

        py::class_<MultipleGridBuilder, LevelGridBuilder, std::shared_ptr<MultipleGridBuilder>>(gridGeneratorModule, "MultipleGridBuilder")
        .def("make_shared", &MultipleGridBuilder::makeShared, py::return_value_policy::reference)
        .def("add_coarse_grid", &MultipleGridBuilder::addCoarseGrid)
        .def("add_grid", py::overload_cast<Object*>(&MultipleGridBuilder::addGrid))
        .def("add_grid", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGrid))
        .def("add_geometry", py::overload_cast<Object*>(&MultipleGridBuilder::addGeometry))
        .def("add_geometry", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGeometry))
        .def("get_number_of_levels", &MultipleGridBuilder::getNumberOfLevels)
        .def("build_grids", &MultipleGridBuilder::buildGrids)
        .def("set_subdomain_box", &MultipleGridBuilder::setSubDomainBox)
        .def("find_communication_indices", &MultipleGridBuilder::findCommunicationIndices)
        .def("set_communication_process", &MultipleGridBuilder::setCommunicationProcess);

        return gridGeneratorModule;
    }
}
