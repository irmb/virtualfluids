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
#ifndef PYBIND_LOGGER
#define PYBIND_LOGGER

#include <pybind11/pybind11.h>
#include <logger/Logger.h>

namespace logger_bindings
{
    namespace py = pybind11;

    inline void makeModule(py::module_ &pyfluids)
    {
        auto logger = pyfluids.def_submodule("logger");
        py::class_<vf::logging::Logger>(logger, "Logger")
        .def_static("initialize_logger", &vf::logging::Logger::initializeLogger)
        .def_static("change_log_path", &vf::logging::Logger::changeLogPath, py::arg("path"));

        // use f-strings (f"text {float}") in python for compounded messages
        logger.def("vf_log_trace", [](const std::string& message){ VF_LOG_TRACE(message); }, py::arg("message"));        
        logger.def("vf_log_debug", [](const std::string& message){ VF_LOG_DEBUG(message); }, py::arg("message"));        
        logger.def("vf_log_info", [](const std::string& message){ VF_LOG_INFO(message); }, py::arg("message"));        
        logger.def("vf_log_warning", [](const std::string& message){ VF_LOG_WARNING(message); }, py::arg("message"));        
        logger.def("vf_log_critical", [](const std::string& message){ VF_LOG_CRITICAL(message); }, py::arg("message"));       
    }
}

#endif
