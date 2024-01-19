//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include "basics/config/ConfigurationFile.h"

namespace configuration
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<vf::basics::ConfigurationFile>(parentModule, "ConfigurationFile")
        .def(py::init<>())
        .def("load", &vf::basics::ConfigurationFile::load, py::arg("file"))
        .def("contains", &vf::basics::ConfigurationFile::contains, py::arg("key"))
        .def("get_int_value"   , static_cast<int         (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"))
        .def("get_int_value"   , static_cast<int         (vf::basics::ConfigurationFile::*)(const std::string&, int        ) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"))
        .def("get_uint_value"  , static_cast<uint        (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"))
        .def("get_uint_value"  , static_cast<uint        (vf::basics::ConfigurationFile::*)(const std::string&, uint       ) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"))
        .def("get_float_value" , static_cast<float       (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"))
        .def("get_float_value" , static_cast<float       (vf::basics::ConfigurationFile::*)(const std::string&, float      ) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"))
        .def("get_double_value", static_cast<double      (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"))
        .def("get_double_value", static_cast<double      (vf::basics::ConfigurationFile::*)(const std::string&, double     ) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"))
        .def("get_bool_value"  , static_cast<bool        (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"))
        .def("get_bool_value"  , static_cast<bool        (vf::basics::ConfigurationFile::*)(const std::string&, bool       ) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"))
        .def("get_string_value", static_cast<std::string (vf::basics::ConfigurationFile::*)(const std::string&) const>(&vf::basics::ConfigurationFile::getValue<std::string>), py::arg("key"))
        .def("get_string_value", static_cast<std::string (vf::basics::ConfigurationFile::*)(const std::string&, std::string) const>(&vf::basics::ConfigurationFile::getValue), py::arg("key"), py::arg("default_value"));
    }
}