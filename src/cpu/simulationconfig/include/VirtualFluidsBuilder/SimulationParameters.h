//
// Created by Sven Marcus on 11.09.20.
//

#ifndef VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
#define VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H

#include <array>
#include <geometry3d/GbPoint3D.h>

struct PhysicalParameters {
    double latticeViscosity{};
    double bulkViscosityFactor{1};
    double latticeDensity{};
};

struct GridParameters {
    std::array<int, 3> numberOfNodesPerDirection{1, 1, 1};
    std::array<int, 3> blocksPerDirection{1, 1, 1};
    int referenceDirectionIndex{};
    double deltaX{1};
    bool periodicBoundaryInX1{};
    bool periodicBoundaryInX2{};
    bool periodicBoundaryInX3{};
};

struct SimulationParameters {
    int numberOfTimeSteps{};
    int timeStepLogInterval{};
    int numberOfThreads{};
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
