#include <pybind11/pybind11.h>
#include "basics/Core/LbmOrGks.h"

namespace lbmOrGks
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
         py::enum_<LbmOrGks>(parentModule, "LbmOrGks")
         .value("LBM", LbmOrGks::LBM)
         .value("GKS", LbmOrGks::GKS);

    }
}