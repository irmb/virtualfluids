#include <pybind11/pybind11.h>
#include "basics/config/ConfigurationFile.h"

namespace configuration
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<vf::basics::ConfigurationFile>(parentModule, "ConfigurationFile")
        .def(py::init<>())
        .def("load", &vf::basics::ConfigurationFile::load)
        .def("contains", &vf::basics::ConfigurationFile::contains)
        .def("get_int_value"   , static_cast<int         (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_int_value"   , static_cast<int         (vf::basics::ConfigurationFile::*)(const std::string&, int        ) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_uint_value"  , static_cast<uint        (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_uint_value"  , static_cast<uint        (vf::basics::ConfigurationFile::*)(const std::string&, uint       ) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_float_value" , static_cast<float       (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_float_value" , static_cast<float       (vf::basics::ConfigurationFile::*)(const std::string&, float      ) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_double_value", static_cast<double      (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_double_value", static_cast<double      (vf::basics::ConfigurationFile::*)(const std::string&, double     ) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_bool_value"  , static_cast<bool        (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_bool_value"  , static_cast<bool        (vf::basics::ConfigurationFile::*)(const std::string&, bool       ) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_string_value", static_cast<std::string (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue))
        .def("get_string_value", static_cast<std::string (vf::basics::ConfigurationFile::*)(const std::string&, std::string) const>(&vf::basics::ConfigurationFile::getValue));
    }
}