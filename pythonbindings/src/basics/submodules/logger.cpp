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
        .def_static("enable_printed_rank_numbers", &logging::Logger::enablePrintedRankNumbers);

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