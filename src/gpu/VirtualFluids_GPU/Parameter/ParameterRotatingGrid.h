#ifndef PARAMETER_ROTATING_GRID
#define PARAMETER_ROTATING_GRID

#include <array>
#include <vector>

#include <basics/DataTypes.h>
#include <basics/geometry3d/Axis.h>

struct ParameterRotatingGrid {
public:
    ParameterRotatingGrid(const std::array<real, 3> &centerPoint, const Axis &rotationalAxis);
    void initializeNestedCoordinates(const std::array<real *, 3> &globalCoordinates, uint numberOfNodes);
    void transformNestedToBase(const std::array<real *, 3> &globalCoordinates);

public:
    const std::array<real, 3> centerPoint;
    const Axis rotationalAxis;

    std::vector<real> nestedCoordinatesX;
    std::vector<real> nestedCoordinatesY;
    std::vector<real> nestedCoordinatesZ;

    real gridAngle = 1.0;
    std::array<real, 3> angularVelocity;
};

#endif
