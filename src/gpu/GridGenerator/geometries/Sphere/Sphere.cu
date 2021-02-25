#include "Sphere.h"

#include <algorithm>    // std::min

#include "geometries/Vertex/Vertex.h"

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

CUDA_HOST int Sphere::getIntersection(const Vertex & point, const Vertex & direction, Vertex & pointOnObject, real & qVal)
{
    
    Vertex relativePoint( point.x - this->centerX, 
                          point.y - this->centerY, 
                          point.z - this->centerZ );

    real directionSquare = direction.x * direction.x
                         + direction.y * direction.y
                         + direction.z * direction.z;

    real p = 2* ( relativePoint.x * direction.x
                + relativePoint.y * direction.y
                + relativePoint.z * direction.z )
                / directionSquare;

    real q = ( relativePoint.x * relativePoint.x
             + relativePoint.y * relativePoint.y
             + relativePoint.z * relativePoint.z
             - (real)this->radius * (real)this->radius )
           / directionSquare;

    real discriminant = 0.25 * p * p - q;


    if( discriminant < 0.0 ) return 1;

    real result1 = - 0.5 * p + std::sqrt( discriminant );
    real result2 = - 0.5 * p - std::sqrt( discriminant );

    if( result1 < 0.0 && result2 < 0.0 ) return 1;

    if (result1 < 0.0)
        result1 = (real)FLT_MAX;
    if (result2 < 0.0)
        result2 = (real)FLT_MAX;

    real t = std::min( result1, result2 );

    pointOnObject.x = point.x + t * direction.x;
    pointOnObject.y = point.y + t * direction.y;
    pointOnObject.z = point.z + t * direction.z;

    qVal = t;

    return 0;
}
