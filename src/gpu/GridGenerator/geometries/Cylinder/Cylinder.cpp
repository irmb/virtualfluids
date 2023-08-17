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

double Cylinder::getRadius() const
{
    return radius;
}

double Cylinder::getHeight() const
{
    return height;
}

Cylinder::PrincipalAxis Cylinder::getPrincipalAxis() const
{
    return principalAxis;
}

bool Cylinder::isInCircle(double delta1, double delta2, double offset) const
{
    return (delta1 * delta1 + delta2 * delta2) < ((this->radius - offset) * (this->radius - offset));
}

bool Cylinder::isPointInObject(const double &x1, const double &x2, const double &x3, const double &minOffset,
                               const double &maxOffset)
{
    double offset = maxOffset;
    if (x1 < center.at(x) || x2 < center.at(y) || x3 < center.at(z)) offset = minOffset;

    const double deltaX1 = x1 - center.at(x);
    const double deltaX2 = x2 - center.at(y);
    const double deltaX3 = x3 - center.at(z);

    switch (principalAxis) {
        case x:
            if (deltaX1 > 0.5 * height || deltaX1 < -0.5 * height) return false;
            return isInCircle(deltaX2, deltaX3, offset);
        case y:
            if (deltaX2 > 0.5 * height || deltaX2 < -0.5 * height) return false;
            return isInCircle(deltaX1, deltaX3, offset);
        case z:
            if (deltaX3 > 0.5 * height || deltaX3 < -0.5 * height) return false;
            return isInCircle(deltaX1, deltaX2, offset);
    }

    VF_LOG_CRITICAL("Unknown principal axis in Cylinder.");
    return false;
}

void Cylinder::changeSizeByDelta(double delta)
{
    this->radius += delta;
    this->height += 2 * delta;
}