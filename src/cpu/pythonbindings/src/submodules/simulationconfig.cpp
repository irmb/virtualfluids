#include <pybind11/pybind11.h>
#include <simulationconfig/Simulation.h>

namespace py = pybind11;

void makeSimulationModule(py::module_ &parentModule)
{
    using namespace pybind11::literals;

    py::class_<Simulation, std::shared_ptr<Simulation>>(parentModule, "Simulation")
            .def(py::init())
            .def("set_writer", &Simulation::setWriterConfig)
            .def("set_grid_parameters", &Simulation::setGridParameters)
            .def("set_physical_parameters", &Simulation::setPhysicalParameters)
            .def("set_simulation_parameters", &Simulation::setSimulationParameters)
            .def("set_kernel_config", &Simulation::setKernelConfig)
            .def("add_object", &Simulation::addObject)
            .def("add_bc_adapter", &Simulation::addBCAdapter)
            .def("run_simulation", &Simulation::run);


}