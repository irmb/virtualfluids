//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBHALFSPACE3D_H
#define GBHALFSPACE3D_H

#include <sstream>
#include <iostream>

#include <basics/utilities/UbMath.h>

#include <numerics/geometry3d/GbPoint3D.h>
#include <numerics/geometry3d/GbTriangle3D.h>
#include <numerics/geometry3d/GbVector3D.h>

#include <basics/memory/MbSharedPointerDefines.h>
class GbHalfSpace3D;
typedef VFSharedPtr<GbHalfSpace3D> GbHalfSpace3DPtr;


/*=========================================================================*/
/* GbHalfSpace3D                                                             */
/*                                                                         */
/**
* This Class helps in performing some operations on a halfspace defined by 2 or 3 points
*/

class GbHalfSpace3D                            
{
public:
   GbHalfSpace3D(GbTriangle3D* triangle);

   GbHalfSpace3D(GbPoint3D* PointA, GbPoint3D* PointB, GbPoint3D* PointC);

   GbHalfSpace3D(GbPoint3D* PointA, GbPoint3D* PointB);

   GbHalfSpace3D(  const double& p1x, const double& p1y, const double& p1z
                 , const double& p2x, const double& p2y, const double& p2z
                 , const double& p3x, const double& p3y, const double& p3z );
   GbHalfSpace3D( const double& p1x, const double& p1y, const double& p1z,
                  const double& nx, const double& ny, const double& nz);

   /*=======================================================*/
   std::string getTypeID() {return "GbHalfSpace3D"; }
   /*=============================================*/
   bool ptInside(const double& x, const double& y, const double& z)
   {
      return UbMath::greaterEqual( normalX*x + normalY*y + normalZ*z, this->d );
   }
   /*=============================================*/
   bool ptInside(GbPoint3D* pointX)
   {
      //GbVector3D X(PointX->x1, PointX->x2, PointX->x3 );
      //return UbMath::greaterEqual(this->Normal.Dot(X), this->d);      
      return UbMath::greaterEqual(  normalX*pointX->x1 + normalY*pointX->x2 + normalZ*pointX->x3, this->d );
   }
   /*=============================================*/
   bool ptInside(GbVector3D& x)
   {
      //return UbMath::greaterEqual(this->Normal.Dot(X), this->d);
      return UbMath::greaterEqual(  normalX*x[0] + normalY*x[1] + normalZ*x[2], this->d );
   }
   /*=============================================*/
   double getDistance(const double& x1p, const double& x2p, const double& x3p)
   {
      return (normalX*x1p + normalY*x2p + normalZ*x3p) - this->d;
   }

   const double& getNormalX() { return this->normalX; }
   const double& getNormalY() { return this->normalY; }
   const double& getNormalZ() { return this->normalZ; }
   const double& getD()       { return this->d;       }

private:
   //GbVector3D Normal;
   double normalX;
   double normalY;
   double normalZ;
   double d;
};
/*=========================================================================*/

#endif //GBHALFSPACE3D_H
