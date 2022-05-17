#include <pybind11/pybind11.h>
#include "gpu/GridGenerator/grid/GridFactory.h"

namespace grid_factory
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {        
        py::class_<GridFactory, std::shared_ptr<GridFactory>>(parentModule, "GridFactory")
        .def("make", &GridFactory::make, py::return_value_policy::reference);
    }
}