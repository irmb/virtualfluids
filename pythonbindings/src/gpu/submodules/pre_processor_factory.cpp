#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactory.h>
#include <gpu/VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactoryImp.h>

namespace pre_processor_factory
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<PreProcessorFactory, std::shared_ptr<PreProcessorFactory>>(parentModule, "_PreProcessorFactory");

        py::class_<PreProcessorFactoryImp, PreProcessorFactory, std::shared_ptr<PreProcessorFactoryImp>>(parentModule, "PreProcessorFactory")
        .def("get_instance", &PreProcessorFactoryImp::getInstance, py::return_value_policy::reference);
    }
}