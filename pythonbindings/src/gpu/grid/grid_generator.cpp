#include <pybind11/pybind11.h>
#include "gpu/GridGenerator/geometries/Object.h"
#include "gpu/GridGenerator/geometries/BoundingBox/BoundingBox.h"
#include "gpu/GridGenerator/geometries/Conglomerate/Conglomerate.h"
#include "gpu/GridGenerator/geometries/Cuboid/Cuboid.h"
#include "gpu/GridGenerator/geometries/Sphere/Sphere.h"
#include "gpu/GridGenerator/geometries/TriangularMesh/TriangularMesh.h"

namespace grid_generator
{
    namespace py = pybind11;
    py::module makeModule(py::module_ &parentModule)
    {  
        py::module gridGeneratorModule = parentModule.def_submodule("grid_generator");

        py::class_<BoundingBox>(gridGeneratorModule, "BoundingBox")
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

        return gridGeneratorModule;
    }
}
