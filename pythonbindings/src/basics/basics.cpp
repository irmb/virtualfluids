#include <pybind11/pybind11.h>
#include "submodules/logger.cpp"
#include "submodules/configuration_file.cpp"

namespace basics
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module basicsModule = parentModule.def_submodule("basics");

        logging::makeModule(basicsModule);
        configuration::makeModule(basicsModule);
        
        return basicsModule;
    }
}