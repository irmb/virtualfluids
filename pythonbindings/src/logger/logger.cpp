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
        loggerModule.def("vf_log_trace", [](std::string message){ VF_LOG_TRACE(message); }, py::arg("message"));        
        loggerModule.def("vf_log_debug", [](std::string message){ VF_LOG_DEBUG(message); }, py::arg("message"));        
        loggerModule.def("vf_log_info", [](std::string message){ VF_LOG_INFO(message); }, py::arg("message"));        
        loggerModule.def("vf_log_warning", [](std::string message){ VF_LOG_WARNING(message); }, py::arg("message"));        
        loggerModule.def("vf_log_critical", [](std::string message){ VF_LOG_CRITICAL(message); }, py::arg("message"));        

        return loggerModule;
    }
} // namespace logging
