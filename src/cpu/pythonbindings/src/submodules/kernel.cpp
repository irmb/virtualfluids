#include <memory>
#include <pybind11/pybind11.h>
#include <simulationconfig/KernelFactory.h>
#include <simulationconfig/KernelConfigStructs.h>


namespace py = pybind11;


void makeKernelModule(py::module &parentModule)
{
    using namespace pybind11::literals;

    py::module kernelModule = parentModule.def_submodule("kernel");

    py::enum_<KernelFactory::KernelType>(kernelModule, "KernelType")
            .value("BGK", KernelFactory::BGK)
            .value("CompressibleCumulantFourthOrderViscosity",
                   KernelFactory::COMPRESSIBLE_CUMULANT_4TH_ORDER_VISCOSITY);


    py::class_<LBMKernelConfig, std::shared_ptr<LBMKernelConfig>>(kernelModule, "LBMKernel")
            .def(py::init<KernelFactory::KernelType>())
            .def_readwrite("use_forcing", &LBMKernelConfig::useForcing)
            .def_readwrite("forcing_in_x1", &LBMKernelConfig::forcingX1)
            .def_readwrite("forcing_in_x2", &LBMKernelConfig::forcingX2)
            .def_readwrite("forcing_in_x3", &LBMKernelConfig::forcingX3)
            .def("set_forcing", [](LBMKernelConfig &kernelConfig, double x1, double x2, double x3) {
                kernelConfig.forcingX1 = x1;
                kernelConfig.forcingX2 = x2;
                kernelConfig.forcingX3 = x3;
            })
            .def("__repr__", [](LBMKernelConfig &kernelConfig) {
                std::ostringstream stream;
                stream << "<" << kernelConfig.kernelType << std::endl
                       << "Use forcing: " << kernelConfig.useForcing << std::endl
                       << "Forcing in x1: " << kernelConfig.forcingX1 << std::endl
                       << "Forcing in x2: " << kernelConfig.forcingX2 << std::endl
                       << "Forcing in x3: " << kernelConfig.forcingX3 << ">" << std::endl;

                return stream.str();
            });
}