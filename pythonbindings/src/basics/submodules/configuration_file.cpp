#include <pybind11/pybind11.h>
#include <basics/config/ConfigurationFile.h>

namespace configuration
{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module configModule = parentModule.def_submodule("configuration");

        py::class_<vf::basics::ConfigurationFile>(configModule, "ConfigurationFile")
        .def("load", &vf::basics::ConfigurationFile::load);

        return configModule;
    }
}