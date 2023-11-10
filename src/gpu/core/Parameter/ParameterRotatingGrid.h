#ifndef PARAMETER_ROTATING_GRID
#define PARAMETER_ROTATING_GRID

#include "PointerDefinitions.h"
#include <array>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/geometry3d/Axis.h>

struct ParameterRotatingGridHostDevice
{
    std::array<real, 3> centerPoint;
    real* nestedCoordinatesX = nullptr; // in local coordinate system of rotating grid
    real* nestedCoordinatesY = nullptr;
    real* nestedCoordinatesZ = nullptr;
    uint sizeOfNestedCoordinates;
    uint memorySizeOfNestedCoordinates;
    std::array<real, 3> gridAngle = { 0.0, 0.0, 0.0 };
    std::array<real, 3> angularVelocity = { 0.0, 0.0, 0.0 };
};

class ParameterRotatingGrid
{
public:
    ParameterRotatingGrid(const std::array<real, 3>& centerPoint, const Axis& rotationalAxis, uint sizeOfNestedCoordinates);
    void fillNestedCoordinateVectorsOnHost(const std::array<real*, 3>& globalCoordinates);

    const Axis rotationalAxis;
    SPtr<ParameterRotatingGridHostDevice> parameterRotHost = std::make_shared<ParameterRotatingGridHostDevice>();
    SPtr<ParameterRotatingGridHostDevice> parameterRotDevice = std::make_shared<ParameterRotatingGridHostDevice>();

    real initialGridRotation = 0.0; // for debugging purposes
};

#endif
