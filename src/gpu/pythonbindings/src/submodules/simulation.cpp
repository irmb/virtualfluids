#include <pybind11/pybind11.h>
#include <VirtualFluids_GPU/LBM/Simulation.h>

namespace simulation{

    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {

        py::module simModule = parentModule.def_submodule("simulation");

        py::class_<Simulation>(simModule, "Simulation")
        .def("set_factories", &Simulation::setFactories)
        .def("init", &Simulation::init)
        .def("run", &Simulation::run)
        .def("free", &Simulation::free);
    }
}