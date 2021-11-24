#include <pybind11/pybind11.h>
#include "submodules/actuator_line.cpp"

namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(pyfluids, m)
    {
        actuator_line::makeModule(m);
    }
}