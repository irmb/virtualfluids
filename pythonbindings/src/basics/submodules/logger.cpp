#include <pybind11/pybind11.h>
#include <basics/Core/Logger/Logger.h>

namespace logging
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module loggerModule = parentModule.def_submodule("logger");

        py::class_<logging::Logger>(loggerModule, "Logger")
        .def("add_stream", &logging::Logger::addStream)
        .def("set_debug_level", &logging::Logger::setDebugLevel)
        .def("time_stamp", &logging::Logger::timeStamp)
        .def("enable_printed_rank_numbers", &logging::Logger::enablePrintedRankNumbers);

        py::enum_<logging::Logger::Level>(loggerModule, "level")
        .value("INFO_LOW", logging::Logger::Level::INFO_LOW)
        .value("INFO_INTERMEDIATE", logging::Logger::Level::INFO_INTERMEDIATE)
        .value("INFO_HIGH", logging::Logger::Level::INFO_HIGH)
        .value("WARNING", logging::Logger::Level::WARNING)
        .value("LOGGER_ERROR", logging::Logger::Level::LOGGER_ERROR);

        py::enum_<logging::Logger::TimeStamp>(loggerModule, "time_stamp")
        .value("ENABLE", logging::Logger::TimeStamp::ENABLE)
        .value("DISABLE", logging::Logger::TimeStamp::DISABLE);

        return loggerModule;
    }
}