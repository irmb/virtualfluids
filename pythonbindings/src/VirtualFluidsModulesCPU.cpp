#include <pybind11/pybind11.h>
#include "cpu/cpu.cpp"

namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(bindings, m)
    {
        cpu::makeModule(m);
    }
}