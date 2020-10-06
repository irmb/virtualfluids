#include <pybind11/pybind11.h>
#include <VirtualFluidsBuilder/VirtualFluidsBuilder.h>

namespace py = pybind11;

void makeBuilderModule(py::module &parentModule)
{
    using namespace pybind11::literals;

    py::module builderModule = parentModule.def_submodule("builder");

    py::class_<VirtualFluidsBuilder, std::shared_ptr<VirtualFluidsBuilder>>(builderModule, "VirtualFluidsBuilder")
            .def(py::init())
            .def("set_writer", &VirtualFluidsBuilder::setWriterConfig)
            .def("set_grid_parameters", &VirtualFluidsBuilder::setGridParameters)
            .def("set_physical_parameters", &VirtualFluidsBuilder::setPhysicalParameters)
            .def("set_simulation_parameters", &VirtualFluidsBuilder::setSimulationParameters)
            .def("set_kernel_config", &VirtualFluidsBuilder::setKernelConfig)
            .def("add_object", &VirtualFluidsBuilder::addObject)
            .def("add_bc_adapter", &VirtualFluidsBuilder::addBCAdapter)
            .def("run_simulation", &VirtualFluidsBuilder::run);


}