#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/GPU/CudaMemoryManager.h>
#include <gpu/VirtualFluids_GPU/Parameter/Parameter.h>


namespace cuda_memory_manager
{

    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<CudaMemoryManager, std::shared_ptr<CudaMemoryManager>>(parentModule, "CudaMemoryManager")
        .def("make", &CudaMemoryManager::make, py::return_value_policy::reference);

    }
}