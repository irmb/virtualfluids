#include <pybind11/pybind11.h>
#include "basics/basics.cpp"
#include "lbm/lbm.cpp"
#include "gpu/gpu.cpp"
#include "logger/logger.cpp"

namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(pyfluids, m)
    {
        basics::makeModule(m);
        gpu::makeModule(m);
        lbm::makeModule(m);
        logger::makeModule(m);
        py::add_ostream_redirect(m, "ostream_redirect");
    }
}