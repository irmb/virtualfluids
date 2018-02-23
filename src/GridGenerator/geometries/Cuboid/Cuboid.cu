#include "Cuboid.h"
#include "utilities/math/CudaMath.cuh"


Cuboid::Cuboid(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b)
    : minX1(x1a), minX2(x2a), minX3(x3a), maxX1(x1b), maxX2(x2b), maxX3(x3b)
{
    x1min = getX1Minimum();
    x2min = getX2Minimum();
    x3min = getX3Minimum();

    x1max = getX1Maximum();
    x2max = getX2Maximum();
    x3max = getX3Maximum();
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
    if (CudaMath::lessEqual(x1, this->getX1Minimum() + minOffset))    return false;
    if (CudaMath::lessEqual(x2, this->getX2Minimum() + minOffset))    return false;
    if (CudaMath::lessEqual(x3, this->getX3Minimum() + minOffset))    return false;
    if (CudaMath::greaterEqual(x1, this->getX1Maximum() - maxOffset)) return false;
    if (CudaMath::greaterEqual(x2, this->getX2Maximum() - maxOffset)) return false;
    if (CudaMath::greaterEqual(x3, this->getX3Maximum() - maxOffset)) return false;

    return true;
}

bool Cuboid::isOnBoundary(const double& x, const double& y, const double& z, const double& minOffset, const double& maxOffset)
{
    const bool isOnXYPlanes = isOn(z, this->minX3 + minOffset, this->maxX3 + maxOffset) && isBetween(y, this->minX2 + minOffset, this->maxX2 + maxOffset) && isBetween(x, this->minX1 + minOffset, this->maxX1 + maxOffset);
    const bool isOnXZPlanes = isOn(y, this->minX2 + minOffset, this->maxX2 + maxOffset) && isBetween(x, this->minX1 + minOffset, this->maxX1 + maxOffset) && isBetween(z, this->minX3 + minOffset, this->maxX3 + maxOffset);
    const bool isOnYZPlanes = isOn(x, this->minX1 + minOffset, this->maxX1 + maxOffset) && isBetween(y, this->minX2 + minOffset, this->maxX2 + maxOffset) && isBetween(z, this->minX3 + minOffset, this->maxX3 + maxOffset);

    return isOnXYPlanes || isOnXZPlanes || isOnYZPlanes;
}


bool Cuboid::isOn(const real& coord, const real& plane1, const real& plane2)
{
    return  CudaMath::equal(coord, plane1) || CudaMath::equal(coord, plane2);
}

bool Cuboid::isBetween(const real& coord, const real& start, const real& end)
{
    return  CudaMath::greaterEqual(coord, start) && CudaMath::lessEqual(coord, end);
}



///*=======================================================*/
//double Cuboid::getLengthX1()
//{ 
//   return (this->getX1Maximum() - this->getX1Minimum() ); 
//}
///*=======================================================*/
//double Cuboid::getLengthX2()
//{ 
//   return (this->getX2Maximum() - this->getX2Minimum());  
//}
///*=======================================================*/
//double Cuboid::getLengthX3()
//{ 
//   return (this->getX3Maximum() - this->getX3Minimum()); 
//}
///*=======================================================*/
//bool Cuboid::isPointInObject(const double& x1p, const double& x2p, const double& x3p)
//{
//   //true, wenn 'in Object' oder 'auf Boundary'!
//   if     (UbMath::less(x1p,this->getX1Minimum()))    return false;
//   else if(UbMath::less(x2p,this->getX2Minimum()))    return false;
//   else if(UbMath::less(x3p,this->getX3Minimum()))    return false;
//   else if(UbMath::greater(x1p,this->getX1Maximum())) return false;
//   else if(UbMath::greater(x2p,this->getX2Maximum())) return false;
//   else if(UbMath::greater(x3p,this->getX3Maximum())) return false;
//
//   return true;
//}
///*=======================================================*/
//bool Cuboid::isPointInObject(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary)
//{
//   pointIsOnBoundary = false;
//   
//   //true, wenn 'in Object' oder 'auf Boundary'!
//   if     (UbMath::less(x1p,this->getX1Minimum()))    return false;
//	else if(UbMath::less(x2p,this->getX2Minimum()))    return false;
//	else if(UbMath::less(x3p,this->getX3Minimum()))    return false;
//	else if(UbMath::greater(x1p,this->getX1Maximum())) return false;
//	else if(UbMath::greater(x2p,this->getX2Maximum())) return false;
//	else if(UbMath::greater(x3p,this->getX3Maximum())) return false;
//	
//   if     (UbMath::equal(x1p,this->getX1Minimum())) pointIsOnBoundary = true;
//   else if(UbMath::equal(x2p,this->getX2Minimum())) pointIsOnBoundary = true;
//   else if(UbMath::equal(x3p,this->getX3Minimum())) pointIsOnBoundary = true;
//   else if(UbMath::equal(x1p,this->getX1Maximum())) pointIsOnBoundary = true;
//   else if(UbMath::equal(x2p,this->getX2Maximum())) pointIsOnBoundary = true;
//   else if(UbMath::equal(x3p,this->getX3Maximum())) pointIsOnBoundary = true;
//   
//   return true;
//}

///*=======================================================*/
//bool Cuboid::isCellInsideGbObject3D(const double& x1p1,const double& x2p1,const double& x3p1,const double& x1p2,const double& x2p2,const double& x3p2)
//{
//   if     ( UbMath::less   (x1p1, this->getX1Minimum() ) ) return false;
//   else if( UbMath::less   (x2p1, this->getX2Minimum() ) ) return false;
//   else if( UbMath::less   (x3p1, this->getX3Minimum() ) ) return false;
//   else if( UbMath::greater(x1p2, this->getX1Maximum() ) ) return false;
//   else if( UbMath::greater(x2p2, this->getX2Maximum() ) ) return false;
//   else if( UbMath::greater(x3p2, this->getX3Maximum() ) ) return false;
//
//   return true;
//}
//
//void Cuboid::translate(const double& tx1, const double& tx2, const double& tx3)
//{  
//   this->p1.translate(tx1, tx2, tx3);
//   this->p2.translate(tx1, tx2, tx3);
//}
///*=======================================================*/
//void Cuboid::scale(const double& sx1, const double& sx2, const double& sx3)
//{  
//   double lenX1 = this->getLengthX1();
//   double lenX2 = this->getLengthX2();
//   double lenX3 = this->getLengthX3();
//
//   double deltaX1 = lenX1*sx1 - lenX1;
//   double deltaX2 = lenX2*sx2 - lenX2;
//   double deltaX3 = lenX3*sx3 - lenX3;
//
//   double p1X1 = this->p1.getX1Coordinate();
//   double p1X2 = this->p1.getX2Coordinate();
//   double p1X3 = this->p1.getX3Coordinate();
//   
//   double p2X1 = this->p2.getX1Coordinate();
//   double p2X2 = this->p2.getX2Coordinate();
//   double p2X3 = this->p2.getX3Coordinate();
//
//   this->p1.setCoordinates(p1X1 - 0.5*deltaX1
//                           ,p1X2 - 0.5*deltaX2
//                           ,p1X3 - 0.5*deltaX3);
//
//   this->p2.setCoordinates(p2X1 + 0.5*deltaX1
//                           ,p2X2 + 0.5*deltaX2
//                           ,p2X3 + 0.5*deltaX3);
//}
