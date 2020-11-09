#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <simulationconfig/SimulationParameters.h>

namespace py = pybind11;

void makeParametersModule(py::module_ &parentModule)
{
    py::module parametersModule = parentModule.def_submodule("parameters");

    py::class_<PhysicalParameters, std::shared_ptr<PhysicalParameters>>(parametersModule, "PhysicalParameters")
            .def(py::init())
            .def_readwrite("bulk_viscosity_factor", &PhysicalParameters::bulkViscosityFactor,
                           "The viscosity of the fluid will be multiplied with this factor to calculate its bulk viscosity. Default is 1.0")
            .def_readwrite("lattice_viscosity", &PhysicalParameters::latticeViscosity, "Lattice viscosity")
            .def_readwrite("lattice_density", &PhysicalParameters::latticeDensity, "Lattice Density");

    py::class_<GridParameters, std::shared_ptr<GridParameters>>(parametersModule, "GridParameters")
            .def(py::init())
            .def_readwrite("delta_x", &GridParameters::deltaX)
            .def_readwrite("reference_direction_index", &GridParameters::referenceDirectionIndex)
            .def_readwrite("number_of_nodes_per_direction", &GridParameters::numberOfNodesPerDirection)
            .def_readwrite("blocks_per_direction", &GridParameters::blocksPerDirection)
            .def_readwrite("periodic_boundary_in_x1", &GridParameters::periodicBoundaryInX1)
            .def_readwrite("periodic_boundary_in_x2", &GridParameters::periodicBoundaryInX2)
            .def_readwrite("periodic_boundary_in_x3", &GridParameters::periodicBoundaryInX3);

    py::class_<SimulationParameters, std::shared_ptr<SimulationParameters>>(parametersModule, "SimulationParameters")
            .def(py::init())
            .def_readwrite("number_of_timesteps", &SimulationParameters::numberOfTimeSteps)
            .def_readwrite("timestep_log_interval", &SimulationParameters::timeStepLogInterval)
            .def_readwrite("number_of_threads", &SimulationParameters::numberOfThreads);

}