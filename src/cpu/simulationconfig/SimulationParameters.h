#ifndef VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
#define VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H

#include <array>
#include <memory>
#include <geometry3d/GbPoint3D.h>

struct PhysicalParameters
{
    double latticeViscosity{};
    double bulkViscosityFactor{1};
};

struct BoundingBox
{
    const double minX1;
    const double minX2;
    const double minX3;
    const double maxX1;
    const double maxX2;
    const double maxX3;

    BoundingBox(double minX1, double minX2, double minX3, double maxX1, double maxX2, double maxX3) :
            minX1(minX1),
            minX2(minX2),
            minX3(minX3),
            maxX1(maxX1),
            maxX2(maxX2),
            maxX3(maxX3)
    {}
};

struct GridParameters
{
    std::array<int, 3> numberOfNodesPerDirection{1, 1, 1};
    std::array<int, 3> blocksPerDirection{1, 1, 1};
    int referenceDirectionIndex{};
    double nodeDistance{1};
    bool periodicBoundaryInX1{};
    bool periodicBoundaryInX2{};
    bool periodicBoundaryInX3{};

    std::shared_ptr<BoundingBox> boundingBox()
    {
        return std::make_shared<BoundingBox>(
                0, 0, 0,
                numberOfNodesPerDirection[0] * nodeDistance,
                numberOfNodesPerDirection[1] * nodeDistance,
                numberOfNodesPerDirection[2] * nodeDistance
        );
    }
};

struct RuntimeParameters
{
    int numberOfTimeSteps{};
    int timeStepLogInterval{};
    int numberOfThreads{};
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_SIMULATIONPARAMETERS_H
