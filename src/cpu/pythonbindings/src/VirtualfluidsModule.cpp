#include <pybind11/pybind11.h>
#include "submodules/boundaryconditions.cpp"
#include "submodules/simulationconfig.cpp"
#include "submodules/geometry.cpp"
#include "submodules/kernel.cpp"
#include "submodules/simulationparameters.cpp"
#include "submodules/writer.cpp"

namespace py = pybind11;

PYBIND11_MODULE(pyfluids, m)
{
    makeBoundaryConditionsModule(m);
    makeSimulationModule(m);
    makeGeometryModule(m);
    makeKernelModule(m);
    makeParametersModule(m);
    makeWriterModule(m);
}