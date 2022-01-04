#include <pybind11/pybind11.h>
#include <logger/Logger.h>


namespace logger
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module loggerModule = parentModule.def_submodule("logger");

        py::class_<vf::logging::Logger>(loggerModule, "Logger")
        .def("initialize_logger", &vf::logging::Logger::initalizeLogger)
        .def("change_log_path", &vf::logging::Logger::changeLogPath);
        return loggerModule;
    }
} // namespace logger
