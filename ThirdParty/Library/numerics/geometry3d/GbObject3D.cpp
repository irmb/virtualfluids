#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/creator/GbObject3DCreator.h>
#include <numerics/geometry3d/GbPoint3D.h>
#include <basics/utilities/UbMath.h>                 

using namespace std;

string GbObject3D::getTypeID()
{
      return this->getCreator()->getTypeID();
}
/*======================================================================*/
bool GbObject3D::isPointInGbObject3D(GbPoint3D* p)
{
   return this->isPointInGbObject3D(p->getX1Centroid(),p->getX2Coordinate(),p->getX3Coordinate());
} 
/*======================================================================*/
bool GbObject3D::isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b) 
{

   if(   this->isPointInGbObject3D(x1a, x2a, x3a) 
      && this->isPointInGbObject3D(x1b, x2a, x3a) 
      && this->isPointInGbObject3D(x1b, x2b, x3a)
      && this->isPointInGbObject3D(x1a, x2b, x3a) 
      && this->isPointInGbObject3D(x1a, x2a, x3b)
      && this->isPointInGbObject3D(x1b, x2a, x3b)
      && this->isPointInGbObject3D(x1b, x2b, x3b)
      && this->isPointInGbObject3D(x1a, x2b, x3b))
   {
      return true;
   }

   return false;
}
/*======================================================================*/
bool GbObject3D::isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b) 
{
   if(   this->isPointInGbObject3D(x1a, x2a, x3a)      
      || this->isPointInGbObject3D(x1b, x2a, x3a)      
      || this->isPointInGbObject3D(x1b, x2b, x3a)      
      || this->isPointInGbObject3D(x1a, x2b, x3a)      
      || this->isPointInGbObject3D(x1a, x2a, x3b)      
      || this->isPointInGbObject3D(x1b, x2a, x3b)      
      || this->isPointInGbObject3D(x1b, x2b, x3b)      
      || this->isPointInGbObject3D(x1a, x2b, x3b) )    
   {
      if(   !this->isPointInGbObject3D(x1a, x2a, x3a) 
         || !this->isPointInGbObject3D(x1b, x2a, x3a) 
         || !this->isPointInGbObject3D(x1b, x2b, x3a) 
         || !this->isPointInGbObject3D(x1a, x2b, x3a)
         || !this->isPointInGbObject3D(x1a, x2a, x3b)
         || !this->isPointInGbObject3D(x1b, x2a, x3b)
         || !this->isPointInGbObject3D(x1b, x2b, x3b)
         || !this->isPointInGbObject3D(x1a, x2b, x3b)) return true;
   }
   return false;
}
/*======================================================================*/
bool GbObject3D::isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   if(   this->isPointInGbObject3D(x1a, x2a, x3a)   
      || this->isPointInGbObject3D(x1b, x2a, x3a)   
      || this->isPointInGbObject3D(x1b, x2b, x3a)   
      || this->isPointInGbObject3D(x1a, x2b, x3a)   
      || this->isPointInGbObject3D(x1a, x2a, x3b)   
      || this->isPointInGbObject3D(x1b, x2a, x3b)   
      || this->isPointInGbObject3D(x1b, x2b, x3b)   
      || this->isPointInGbObject3D(x1a, x2b, x3b))  
   {
      return true;
   }

   return false;
}
/*=======================================================*/
bool GbObject3D::isInsideCell(const double& minX1,const double& minX2,const double& minX3,const double& maxX1,const double& maxX2,const double& maxX3)
{
   if(   UbMath::greaterEqual(this->getX1Minimum(),minX1)
      && UbMath::greaterEqual(this->getX2Minimum(),minX2)
      && UbMath::greaterEqual(this->getX3Minimum(),minX3)
      && UbMath::lessEqual(this->getX1Maximum(),maxX1) 
      && UbMath::lessEqual(this->getX2Maximum(),maxX2)   
      && UbMath::lessEqual(this->getX2Maximum(),maxX3)     ) return true;

   return false;
}


