#include "pybind11/pybind11.h"
#include "gpu/VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"
#include "gpu/VirtualFluids_GPU/LBM/LB.h"

namespace turbulence_model
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::enum_<TurbulenceModel>(parentModule, "TurbulenceModel")
        .value("Smagorinsky", TurbulenceModel::Smagorinsky)
        .value("AMD", TurbulenceModel::AMD)
        .value("QR", TurbulenceModel::QR)
        .value("None", TurbulenceModel::None);

        py::class_<TurbulenceModelFactory, std::shared_ptr<TurbulenceModelFactory>>(parentModule, "TurbulenceModelFactory")
        .def(py::init< std::shared_ptr<Parameter>>(), "para")
        .def("set_turbulence_model", &TurbulenceModelFactory::setTurbulenceModel)
        .def("set_model_constant", &TurbulenceModelFactory::setModelConstant)
        .def("read_config_file", &TurbulenceModelFactory::readConfigFile);

    }
}