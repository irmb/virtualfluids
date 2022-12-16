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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file logger.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <pybind11/iostream.h>
#include <basics/Core/Logger/Logger.h>
#include <basics/Core/Logger/implementations/LoggerImp.h>

namespace logger
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module loggerModule = parentModule.def_submodule("logger");

        py::class_<logging::Logger>(loggerModule, "Logger")
        .def_static("add_stdout", [](){
            logging::Logger::addStream(&std::cout);
        })
        .def_static("set_debug_level", &logging::Logger::setDebugLevel)
        .def_static("time_stamp", &logging::Logger::timeStamp, py::arg("time_stamp"))
        .def_static("enable_printed_rank_numbers", &logging::Logger::enablePrintedRankNumbers, py::arg("print"));

        loggerModule.attr("log") = logging::out;
        py::enum_<logging::Logger::Level>(loggerModule, "Level")
        .value("INFO_LOW", logging::Logger::Level::INFO_LOW)
        .value("INFO_INTERMEDIATE", logging::Logger::Level::INFO_INTERMEDIATE)
        .value("INFO_HIGH", logging::Logger::Level::INFO_HIGH)
        .value("WARNING", logging::Logger::Level::WARNING)
        .value("LOGGER_ERROR", logging::Logger::Level::LOGGER_ERROR);

        py::enum_<logging::Logger::TimeStamp>(loggerModule, "TimeStamp")
        .value("ENABLE", logging::Logger::TimeStamp::ENABLE)
        .value("DISABLE", logging::Logger::TimeStamp::DISABLE);

        return loggerModule;
    }
}