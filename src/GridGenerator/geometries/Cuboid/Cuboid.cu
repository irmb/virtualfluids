#include "Cuboid.h"
#include "utilities/math/Math.h"


Cuboid::Cuboid(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b)
    : minX1(x1a), minX2(x2a), minX3(x3a), maxX1(x1b), maxX2(x2b), maxX3(x3b)
{

}

Cuboid::~Cuboid()
{

}

Object* Cuboid::clone() const
{
    return new Cuboid(minX1, minX2, minX3, maxX1, maxX2, maxX3);
}

double Cuboid::getX1Centroid()
{
   return getCenter(minX1, maxX1);
}

double Cuboid::getX1Minimum()
{
    return getMinimum(minX1, maxX1);
}

double Cuboid::getX1Maximum()
{
    return getMaximum(minX1, maxX1);
}

double Cuboid::getX2Centroid()
{
    return getCenter(minX2, maxX2);
}

double Cuboid::getX2Minimum()
{
    return getMinimum(minX2, maxX2);
}	

double Cuboid::getX2Maximum()
{
    return getMaximum(minX2, maxX2);
}

double Cuboid::getX3Centroid()
{
    return getCenter(minX3, maxX3);
}

double Cuboid::getX3Minimum()
{	
    return getMinimum(minX3, maxX3);
}	

double Cuboid::getX3Maximum()
{
    return getMaximum(minX3, maxX3);
}

double Cuboid::getCenter(double x1, double x2)
{
    return 0.5 * (x1 + x2);
}

double Cuboid::getMinimum(double x1, double x2)
{
    return (x1 < x2 ? x1 : x2);
}

double Cuboid::getMaximum(double x1, double x2)
{
    return (x1 > x2 ? x1 : x2);
}

bool Cuboid::isPointInObject(const double& x1, const double& x2, const double& x3, const double& minOffset, const double& maxOffset)
{
    //false, if 'not in Object' or 'on Boundary'!
    if (vf::Math::lessEqual(x1, this->getX1Minimum() + minOffset))    return false;
    if (vf::Math::lessEqual(x2, this->getX2Minimum() + minOffset))    return false;
    if (vf::Math::lessEqual(x3, this->getX3Minimum() + minOffset))    return false;
    if (vf::Math::greaterEqual(x1, this->getX1Maximum() - maxOffset)) return false;
    if (vf::Math::greaterEqual(x2, this->getX2Maximum() - maxOffset)) return false;
    if (vf::Math::greaterEqual(x3, this->getX3Maximum() - maxOffset)) return false;

    return true;
}


bool Cuboid::isOn(const real& coord, const real& plane1, const real& plane2)
{
    return  vf::Math::equal(coord, plane1) || vf::Math::equal(coord, plane2);
}

bool Cuboid::isBetween(const real& coord, const real& start, const real& end)
{
    return  vf::Math::greaterEqual(coord, start) && vf::Math::lessEqual(coord, end);
}


void Cuboid::scale(double delta)
{
    this->minX1 -= delta;
    this->minX2 -= delta;
    this->minX3 -= delta;
                   
    this->maxX1 += delta;
    this->maxX2 += delta;
    this->maxX3 += delta;
}
