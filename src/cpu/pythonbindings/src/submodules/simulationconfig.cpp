#include <pybind11/pybind11.h>
#include <simulationconfig/Simulation.h>

namespace simulation
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<Simulation, std::shared_ptr<Simulation>>(parentModule, "Simulation")
                .def(py::init())
                .def("set_writer", &Simulation::setWriterConfiguration)
                .def("set_grid_parameters", &Simulation::setGridParameters)
                .def("set_physical_parameters", &Simulation::setPhysicalParameters)
                .def("set_runtime_parameters", &Simulation::setRuntimeParameters)
                .def("set_kernel_config", &Simulation::setKernelConfiguration)
                .def("add_object", &Simulation::addObject)
                .def("add_bc_adapter", &Simulation::addBCAdapter)
                .def("run_simulation", &Simulation::run);
    }

}