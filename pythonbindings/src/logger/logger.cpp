#include <pybind11/pybind11.h>
#include <logger/Logger.h>

namespace logging
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module loggerModule = parentModule.def_submodule("logger");

        py::class_<vf::logging::Logger>(loggerModule, "Logger")
        .def("initialize_logger", &vf::logging::Logger::initalizeLogger)
        .def("change_log_path", &vf::logging::Logger::changeLogPath);

        // use f-strings (f"text {float}") in python for compounded messages
        loggerModule.def("vf_log_trace", [](std::string arg){ VF_LOG_TRACE(arg); });        
        loggerModule.def("vf_log_debug", [](std::string arg){ VF_LOG_DEBUG(arg); });        
        loggerModule.def("vf_log_info", [](std::string arg){ VF_LOG_INFO(arg); });        
        loggerModule.def("vf_log_warning", [](std::string arg){ VF_LOG_WARNING(arg); });        
        loggerModule.def("vf_log_critical", [](std::string arg){ VF_LOG_CRITICAL(arg); });        

        return loggerModule;
    }
} // namespace logging
