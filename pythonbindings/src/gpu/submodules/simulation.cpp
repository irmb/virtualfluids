#include <pybind11/pybind11.h>
#include <gpu/VirtualFluids_GPU/LBM/Simulation.h>
#include <gpu/VirtualFluids_GPU/Communication/Communicator.h>
#include <gpu/VirtualFluids_GPU/Kernel/Utilities/KernelFactory/KernelFactory.h>
#include <gpu/VirtualFluids_GPU/PreProcessor/PreProcessorFactory/PreProcessorFactory.h>
#include <gpu/VirtualFluids_GPU/DataStructureInitializer/GridProvider.h>
#include <gpu/VirtualFluids_GPU/Parameter/Parameter.h>
#include <gpu/VirtualFluids_GPU/GPU/CudaMemoryManager.h>
#include <gpu/VirtualFluids_GPU/DataStructureInitializer/GridProvider.h>
#include <gpu/VirtualFluids_GPU/Output/DataWriter.h>

namespace simulation
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<Simulation>(parentModule, "Simulation")
        .def(py::init<>())
        .def("set_factories", &Simulation::setFactories)
        .def("init", &Simulation::init)
        .def("run", &Simulation::run)
        .def("free", &Simulation::free);
    }
}