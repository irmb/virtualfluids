#include <pybind11/pybind11.h>
#include "submodules/boundaryconditions.cpp"
#include "submodules/simulationconfig.cpp"
#include "submodules/geometry.cpp"
#include "submodules/kernel.cpp"
#include "submodules/simulationparameters.cpp"
#include "submodules/writer.cpp"

namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(pyfluids, m)
    {
        boundaryconditions::makeModule(m);
        simulation::makeModule(m);
        geometry::makeModule(m);
        kernel::makeModule(m);
        parameters::makeModule(m);
        writer::makeModule(m);
    }
}