#include <pybind11/pybind11.h>

namespace lbm
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module lbmModule = parentModule.def_submodule("lbm");

        return lbmModule;
    }
}