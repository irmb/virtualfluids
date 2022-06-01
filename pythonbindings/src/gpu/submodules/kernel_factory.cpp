#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactoryImp.h>
#include <gpu/VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactory.h>

namespace kernel_factory
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<KernelFactory, std::shared_ptr<KernelFactory>>(parentModule, "_KernelFactory");
        
        py::class_<KernelFactoryImp, KernelFactory, std::shared_ptr<KernelFactoryImp>>(parentModule, "KernelFactory")
        .def("get_instance", &KernelFactoryImp::getInstance, py::return_value_policy::reference);
    }
}