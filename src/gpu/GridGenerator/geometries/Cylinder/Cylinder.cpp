#include "Cylinder.h"
#include <numeric>

Cylinder::Cylinder(double centerX, double centerY, double centerZ, double radius, double height, PrincipalAxis axis)
    : center({ centerX, centerY, centerZ }), radius(radius), height(height), principalAxis(axis)
{
}

Cylinder::Cylinder(std::array<double, 3> center, double radius, double height, PrincipalAxis axis)
    : center(center), radius(radius), height(height), principalAxis(axis)
{
}

SPtr<Object> Cylinder::clone() const
{
    return std::make_shared<Cylinder>(center, radius, height, principalAxis);
}

double Cylinder::getCentroidCoordinate(PrincipalAxis coordinateDirection) const
{
    return center.at(coordinateDirection);
}

double Cylinder::getMinimunCoordinate(PrincipalAxis coordinateDirection) const
{
    const auto unitVector = unitVectors.at(principalAxis);
    return center.at(coordinateDirection) - 0.5 * height * unitVector.at(coordinateDirection) +
           radius * (unitVector.at(coordinateDirection) - 1);
}

double Cylinder::getMaximumCoordinate(PrincipalAxis coordinateDirection) const
{
    const auto unitVector = unitVectors.at(principalAxis);
    return center.at(coordinateDirection) + 0.5 * height * unitVector.at(coordinateDirection) -
           radius * (unitVector.at(coordinateDirection) - 1);
}

double Cylinder::getX1Centroid()
{
    return getCentroidCoordinate(x);
}

double Cylinder::getX1Minimum()
{
    return getMinimunCoordinate(x);
}

double Cylinder::getX1Maximum()
{
    return getMaximumCoordinate(x);
}

double Cylinder::getX2Centroid()
{
    return getCentroidCoordinate(y);
}

double Cylinder::getX2Minimum()
{
    return getMinimunCoordinate(y);
}

double Cylinder::getX2Maximum()
{
    return getMaximumCoordinate(y);
}

double Cylinder::getX3Centroid()
{
    return getCentroidCoordinate(z);
}

double Cylinder::getX3Minimum()
{
    return getMinimunCoordinate(z);
}

double Cylinder::getX3Maximum()
{
    return getMaximumCoordinate(z);
}

////////////////// #TODO

bool Cylinder::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
 return false;
}


void Cylinder::scale(double delta)
{
}