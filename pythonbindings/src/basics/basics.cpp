#include <pybind11/pybind11.h>
#include "submodules/logger.cpp"
#include "submodules/configuration_file.cpp"
#include "submodules/lbm_or_gks.cpp"

namespace basics
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module basicsModule = parentModule.def_submodule("basics");

        logger::makeModule(basicsModule);
        configuration::makeModule(basicsModule);
        lbmOrGks::makeModule(basicsModule);
        
        return basicsModule;
    }
}