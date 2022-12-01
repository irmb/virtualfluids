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
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"),
                        py::arg("gridScalingFactory"))
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::gpu::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*>(), 
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"))
        .def(py::init<  std::shared_ptr<Parameter>,
                        std::shared_ptr<CudaMemoryManager>,
                        vf::gpu::Communicator &,
                        GridProvider &,
                        BoundaryConditionFactory*,
                        std::shared_ptr<TurbulenceModelFactory>,
                        GridScalingFactory*>(), 
                        py::arg("parameter"),
                        py::arg("memoryManager"),
                        py::arg("communicator"),
                        py::arg("gridProvider"),
                        py::arg("bcFactory"),
                        py::arg("tmFactory"),
                        py::arg("gridScalingFactory"))
        .def("run", &Simulation::run)
        .def("addKineticEnergyAnalyzer", &Simulation::addKineticEnergyAnalyzer, py::arg("t_analyse"))
        .def("addEnstrophyAnalyzer", &Simulation::addEnstrophyAnalyzer, py::arg("t_analyse"));
    }
}