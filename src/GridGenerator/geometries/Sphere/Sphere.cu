#include "Sphere.h"
#include "utilities/math/CudaMath.cuh"

Sphere::Sphere(const double& centerX, const double& centerY, const double& centerZ, const double& radius)
    : centerX(centerX), centerY(centerY), centerZ(centerZ), radius(radius)
{

}

Sphere::~Sphere()
{
}

SPtr<Sphere> Sphere::makeShared(double centerX, double centerY, double centerZ, double radius)
{
    return SPtr<Sphere>(new Sphere(centerX, centerY, centerZ, radius));
}

Object* Sphere::clone() const
{
    return new Sphere(centerX, centerY, centerZ, radius);
}

double Sphere::getX1Centroid()
{
    return centerX;
}

double Sphere::getX1Minimum()
{
    return centerX - radius;
}

double Sphere::getX1Maximum()
{
    return centerX + radius;
}

double Sphere::getX2Centroid()
{
    return centerY;
}

double Sphere::getX2Minimum()
{
    return centerY - radius;
}

double Sphere::getX2Maximum()
{
    return centerY + radius;
}

double Sphere::getX3Centroid()
{
    return centerZ;
}

double Sphere::getX3Minimum()
{
    return centerZ - radius;
}

double Sphere::getX3Maximum()
{
    return centerZ + radius;
}

bool Sphere::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset,
    const double& maxOffset)
{
    double offset = maxOffset;
    if (x1 < centerX || x2 < centerY || x3 < centerZ)
        offset = minOffset;
        

    const double deltaX1 = x1 - centerX;
    const double deltaX2 = x2 - centerY;
    const double deltaX3 = x3 - centerZ;

    return (deltaX1*deltaX1 + deltaX2*deltaX2 + deltaX3*deltaX3) < ((this->radius - offset) * (this->radius - offset));
}


void Sphere::scale(double delta)
{
    this->radius += delta;
}
