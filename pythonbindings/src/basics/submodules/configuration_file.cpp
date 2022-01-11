#include <pybind11/pybind11.h>
#include <basics/config/ConfigurationFile.h>

namespace configuration
{

    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {

        py::class_<vf::basics::ConfigurationFile>(parentModule, "ConfigurationFile")
        .def(py::init<>())
        .def("load", &vf::basics::ConfigurationFile::load);

    }
}