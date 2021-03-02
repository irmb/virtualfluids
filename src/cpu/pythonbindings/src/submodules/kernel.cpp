#include <memory>
#include <pybind11/pybind11.h>
#include <simulationconfig/KernelFactory.h>
#include <simulationconfig/KernelConfigStructs.h>

namespace kernel
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::module kernelModule = parentModule.def_submodule("kernel");

        py::enum_<KernelFactory::KernelType>(kernelModule, "KernelType")
                .value("BGK", KernelFactory::BGK)
                .value("CompressibleCumulantFourthOrderViscosity",
                       KernelFactory::COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY);


        py::class_<LBMKernelConfiguration, std::shared_ptr<LBMKernelConfiguration>>(kernelModule, "LBMKernel")
                .def(py::init<KernelFactory::KernelType>())
                .def_readwrite("use_forcing", &LBMKernelConfiguration::useForcing)
                .def_readwrite("forcing_in_x1", &LBMKernelConfiguration::forcingX1)
                .def_readwrite("forcing_in_x2", &LBMKernelConfiguration::forcingX2)
                .def_readwrite("forcing_in_x3", &LBMKernelConfiguration::forcingX3)
                .def("set_forcing", [](LBMKernelConfiguration &kernelConfig, double x1, double x2, double x3)
                {
                    kernelConfig.forcingX1 = x1;
                    kernelConfig.forcingX2 = x2;
                    kernelConfig.forcingX3 = x3;
                })
                .def("__repr__", [](LBMKernelConfiguration &kernelConfig)
                {
                    std::ostringstream stream;
                    stream << "<" << kernelConfig.kernelType << std::endl
                           << "Use forcing: " << kernelConfig.useForcing << std::endl
                           << "Forcing in x1: " << kernelConfig.forcingX1 << std::endl
                           << "Forcing in x2: " << kernelConfig.forcingX2 << std::endl
                           << "Forcing in x3: " << kernelConfig.forcingX3 << ">" << std::endl;

                    return stream.str();
                });
    }

}