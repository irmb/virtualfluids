#include <pybind11/pybind11.h>
#include "gpu/VirtualFluids_GPU/DataStructureInitializer/GridProvider.h"
// #include <gpu/VirtualFluids_GPU/GPU/CudaMemoryManager.h>
// #include <gpu/VirtualFluids_GPU/Parameter/Parameter.h>
// #include "gpu/GridGenerator/grid/GridBuilder/GridBuilder.h"

namespace grid_provider
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<GridProvider, std::shared_ptr<GridProvider>>(parentModule, "GridProvider")
        .def("make_grid_generator", &GridProvider::makeGridGenerator, py::return_value_policy::reference);
    }
}