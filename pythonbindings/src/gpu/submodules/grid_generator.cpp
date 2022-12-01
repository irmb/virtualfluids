#include <pybind11/pybind11.h>
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

        py::class_<GridFactory, std::shared_ptr<GridFactory>>(gridGeneratorModule, "GridFactory")
        .def_static("make", &GridFactory::make, py::return_value_policy::reference);

        py::class_<BoundingBox>(gridGeneratorModule, "BoundingBox")
        .def(py::init<real, real, real, real, real, real>(), py::arg("min_x"), py::arg("max_x"), py::arg("min_y"), py::arg("max_y"), py::arg("min_z"), py::arg("max_z"));

        py::class_<Object, std::shared_ptr<Object>>(gridGeneratorModule, "Object");
        
        py::class_<Conglomerate, Object, std::shared_ptr<Conglomerate>>(gridGeneratorModule, "Conglomerate")
        .def_static("make_shared", &Conglomerate::makeShared, py::return_value_policy::reference)
        .def("add", &Conglomerate::add, py::arg("object"))
        .def("subtract", &Conglomerate::subtract, py::arg("object"));

        py::class_<Cuboid, Object, std::shared_ptr<Cuboid>>(gridGeneratorModule, "Cuboid")
        .def(py::init<const double&, const double&, const double&, const double&, const double&, const double&>(),
                        py::arg("min_x1"), py::arg("min_x2"), py::arg("min_x3"), py::arg("max_x1"), py::arg("max_x2"), py::arg("max_x3"));

        py::class_<Sphere, Object, std::shared_ptr<Sphere>>(gridGeneratorModule, "Sphere")
        .def_static("make_shared", &Sphere::makeShared, py::return_value_policy::reference);

        py::class_<TriangularMesh, Object, std::shared_ptr<TriangularMesh>>(gridGeneratorModule, "TriangularMesh")
        .def_static("make", &TriangularMesh::make, py::return_value_policy::reference);

        py::class_<GridBuilder, std::shared_ptr<GridBuilder>>(gridGeneratorModule, "GridBuilder")
        .def("get_number_of_grid_levels", &GridBuilder::getNumberOfGridLevels)
        .def("get_grid", &GridBuilder::getGrid, py::arg("level"));

        py::class_<LevelGridBuilder, GridBuilder, std::shared_ptr<LevelGridBuilder>>(gridGeneratorModule, "LevelGridBuilder")
        .def("get_grid", py::overload_cast<int, int>(&LevelGridBuilder::getGrid), py::arg("level"), py::arg("box"))
        .def("set_slip_boundary_condition", &LevelGridBuilder::setSlipBoundaryCondition, py::arg("side_type"), py::arg("normal_x"), py::arg("normal_y"), py::arg("normal_z"))
        .def("set_velocity_boundary_condition", &LevelGridBuilder::setVelocityBoundaryCondition, py::arg("side_type"), py::arg("vx"), py::arg("vy"), py::arg("vz"))
        .def("set_pressure_boundary_condition", &LevelGridBuilder::setPressureBoundaryCondition, py::arg("side_type"), py::arg("rho"))
        .def("set_periodic_boundary_condition", &LevelGridBuilder::setPeriodicBoundaryCondition, py::arg("periodic_x"), py::arg("periodic_y"), py::arg("periodic_z"))
        .def("set_no_slip_boundary_condition", &LevelGridBuilder::setNoSlipBoundaryCondition, py::arg("side_type"))
        .def("set_precursor_boundary_condition", &LevelGridBuilder::setPrecursorBoundaryCondition, py::arg("side_type"), py::arg("file_collection"), py::arg("n_t_read"), py::arg("velocity_x")=0.0f, py::arg("velocity_y")=0.0f, py::arg("velocity_z")=0.0f, py::arg("file_level_to_grid_level_map")=std::vector<uint>())
        .def("set_stress_boundary_condition", &LevelGridBuilder::setStressBoundaryCondition, py::arg("side_type"), py::arg("normal_x"), py::arg("normal_y"), py::arg("normal_z"), py::arg("sampling_offset"), py::arg("z0"), py::arg("dx"));

        py::class_<MultipleGridBuilder, LevelGridBuilder, std::shared_ptr<MultipleGridBuilder>>(gridGeneratorModule, "MultipleGridBuilder")
        .def_static("make_shared", &MultipleGridBuilder::makeShared, py::return_value_policy::reference, py::arg("grid_factory"))
        .def("add_coarse_grid", &MultipleGridBuilder::addCoarseGrid, py::arg("start_x"), py::arg("start_y"), py::arg("start_z"), py::arg("end_x"), py::arg("end_y"), py::arg("end_z"), py::arg("delta"))
        .def("add_grid", py::overload_cast<Object*>(&MultipleGridBuilder::addGrid), py::arg("grid_shape"))
        .def("add_grid", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGrid), py::arg("grid_shape"), py::arg("level_fine"))
        .def("add_geometry", py::overload_cast<Object*>(&MultipleGridBuilder::addGeometry), py::arg("solid_object"))
        .def("add_geometry", py::overload_cast<Object*, uint>(&MultipleGridBuilder::addGeometry), py::arg("solid_object"), py::arg("level"))
        .def("get_number_of_levels", &MultipleGridBuilder::getNumberOfLevels)
        .def("build_grids", &MultipleGridBuilder::buildGrids, py::arg("lbm_or_gks"), py::arg("enable_thin_walls"));

        return gridGeneratorModule;
    }
}
