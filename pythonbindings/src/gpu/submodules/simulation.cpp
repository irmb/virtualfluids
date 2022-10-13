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
#include "gpu/VirtualFluids_GPU/Factories/BoundaryConditionFactory.h"
#include "gpu/VirtualFluids_GPU/TurbulenceModels/TurbulenceModelFactory.h"
#include "gpu/VirtualFluids_GPU/Factories/GridScalingFactory.h"

namespace simulation
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        // missing setFactories and setDataWriter, not possible to wrap these functions as long as they take unique ptr arguments
        py::class_<Simulation>(parentModule, "Simulation")
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::gpu::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*,
                        GridScalingFactory*>(), 
                        "parameter",
                        "memoryManager",
                        "communicator",
                        "gridProvider",
                        "bcFactory",
                        "gridScalingFactory")
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::gpu::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*>(), 
                        "parameter",
                        "memoryManager",
                        "communicator",
                        "gridProvider",
                        "bcFactory")
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::gpu::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*,
                        std::shared_ptr<TurbulenceModelFactory>,
                        GridScalingFactory*>(), 
                        "parameter",
                        "memoryManager",
                        "communicator",
                        "gridProvider",
                        "bcFactory",
                        "tmFactory",
                        "gridScalingFactory")
        .def("run", &Simulation::run)
        .def("addKineticEnergyAnalyzer", &Simulation::addKineticEnergyAnalyzer)
        .def("addEnstrophyAnalyzer", &Simulation::addEnstrophyAnalyzer);
    }
}