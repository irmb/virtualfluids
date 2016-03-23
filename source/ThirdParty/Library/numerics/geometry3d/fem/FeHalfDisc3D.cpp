#include <numerics/geometry3d/fem/FeHalfDisc3D.h>    
#include <numerics/geometry3d/GbSystem3D.h>
#include <numerics/geometry3d/GbLine3D.h>
#include <numerics/geometry3d/GbTriangle3D.h>
#include <numerics/geometry3d/fem/FePoint3D.h>
#include <basics/utilities/UbInfinity.h>

using namespace std;

/*=======================================================*/
ObObjectCreator* FeHalfDisc3D::getCreator()
{
   return NULL;//FeHalfDisc3DCreator::getInstance(); 
}
// Konstruktor
/*==========================================================*/
FeHalfDisc3D::FeHalfDisc3D()
{
   GbPoint3D* p1 = new GbPoint3D();
   GbPoint3D* p2 = new GbPoint3D();
   mLine = new GbLine3D(p1,p2);
   this->mLine->addObserver(this);
   mRad = 0.0;
   cylinderType = FeHalfDisc3D::NOTPARALLELTOAXIS;
}                                                   
/*=======================================================*/
FeHalfDisc3D::FeHalfDisc3D(FeHalfDisc3D* cylinder)
{
   mRad         = cylinder->getRadius();
   cylinderType = cylinder->cylinderType;
   mLine        = cylinder->getLine()->clone();

   this->mLine->addObserver(this);
}                                                   
/*==========================================================*/
FeHalfDisc3D::FeHalfDisc3D(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b, const double& rad)
{
   mLine = new GbLine3D;
   mLine->setPoints( new GbPoint3D(min(x1a,x1b), min(x2a,x2b), min(x3a,x3b))
	                 ,new GbPoint3D(max(x1a,x1b), max(x2a,x2b), max(x3a,x3b)));
   this->mLine->addObserver(this);
   mRad = fabs(rad);

   this->initCylinderType();
   if((this->cylinderType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}                                                            
/*==========================================================*/
FeHalfDisc3D::FeHalfDisc3D(GbPoint3D* p1, GbPoint3D* p2, const double& rad)
{
   mRad = rad;

   mLine = new GbLine3D(p1,p2);
   this->mLine->addObserver(this);
   this->initCylinderType();
   if((this->cylinderType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
FeHalfDisc3D::FeHalfDisc3D(GbLine3D* line, const double& rad)
{
   mRad = rad;

   this->mLine = line;
   this->mLine->addObserver(this);
   this->initCylinderType();
   if((this->cylinderType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/		
// Destruktor
FeHalfDisc3D::~FeHalfDisc3D()
{
   if(mLine) this->mLine->removeObserver(this);
   mLine = NULL;
}
/*=======================================================*/
void FeHalfDisc3D::initCylinderType()
{
   double x1a = mLine->getPoint1()->x1;    double x1b = mLine->getPoint2()->x1;
   double x2a = mLine->getPoint1()->x2;    double x2b = mLine->getPoint2()->x2;
   double x3a = mLine->getPoint1()->x3;    double x3b = mLine->getPoint2()->x3;
   
   if     (x1a!=x1b && x2a==x2b && x3a==x3b)  this->cylinderType = X1PARALLEL; 
   else if(x2a!=x2b && x1a==x1b && x3a==x3b)  this->cylinderType = X2PARALLEL; 
   else if(x3a!=x3b && x1a==x1b && x2a==x2b)  this->cylinderType = X3PARALLEL; 
   else                                       this->cylinderType = NOTPARALLELTOAXIS;
}
/*=======================================================*/
void FeHalfDisc3D::finalize() 
{ 
   if(this->mLine) 
   {
      mLine->finalize();
      delete mLine; 
      mLine=NULL;
   } 
}
/*=======================================================*/
double FeHalfDisc3D::getHeight()
{
   if(mLine) return mLine->getLength(); return 0.0; 
}
/*=======================================================*/
GbPoint3D* FeHalfDisc3D::getPoint1()
{
   if(this->mLine) return this->mLine->getPoint1();
   return NULL;
}
/*=======================================================*/
GbPoint3D* FeHalfDisc3D::getPoint2()
{
   if(this->mLine) return this->mLine->getPoint2();
   return NULL;
}
/*=======================================================*/
void FeHalfDisc3D::setRadius(const double& radius) 
{ 
   this->mRad = std::fabs(radius); 
   this->notifyObserversObjectChanged();
}
/*=======================================================*/
void FeHalfDisc3D::setLine(GbLine3D* line) 
{
   if(this->mLine) this->mLine->removeObserver(this);
   this->mLine = line;  
   this->mLine->addObserver(this);
   this->initCylinderType();

   this->notifyObserversObjectChanged();
}
/*=======================================================*/
void FeHalfDisc3D::setPoint1(const double& x1, const double& x2, const double& x3)
{ 
   if(!mLine->getPoint1()) throw UbException(UB_EXARGS,"line has no point1");
   mLine->getPoint1()->setCoordinates(x1,x2,x3);
   this->initCylinderType();

   //this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*=======================================================*/
void FeHalfDisc3D::setPoint2(const double& x1, const double& x2, const double& x3)
{ 
   if(!mLine->getPoint2()) throw UbException(UB_EXARGS,"line has no point2");
   mLine->getPoint2()->setCoordinates(x1,x2,x3);
   this->initCylinderType();

   //this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*==========================================================*/		
double FeHalfDisc3D::getX1Centroid()  
{
   return mLine->getX1Centroid();
}
/*==========================================================*/
double FeHalfDisc3D::getX1Minimum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX1Minimum(); 
   else if(this->isParallelToX2Axis()) return mLine->getX1Centroid()-mRad; 
   else if(this->isParallelToX3Axis()) return mLine->getX1Centroid()-mRad; 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeHalfDisc3D::getX1Maximum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX1Maximum(); 
   else if(this->isParallelToX2Axis()) return mLine->getX1Centroid()+mRad; 
   else if(this->isParallelToX3Axis()) return mLine->getX1Centroid()+mRad; 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeHalfDisc3D::getX2Centroid()
{
   return mLine->getX2Centroid();
}
/*==========================================================*/
double FeHalfDisc3D::getX2Minimum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX2Centroid()-mRad;
   else if(this->isParallelToX2Axis()) return mLine->getX2Minimum();
   else if(this->isParallelToX3Axis()) return mLine->getX2Centroid()-mRad; 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}	
/*==========================================================*/
double FeHalfDisc3D::getX2Maximum()   
{
   if     (this->isParallelToX1Axis())  return mLine->getX2Centroid()+mRad;
   else if(this->isParallelToX2Axis())  return mLine->getX2Maximum();
   else if(this->isParallelToX3Axis())  return mLine->getX2Centroid()+mRad; 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeHalfDisc3D::getX3Centroid()
{
   return mLine->getX3Centroid();
}
/*==========================================================*/
double FeHalfDisc3D::getX3Minimum()   
{	
   if     (this->isParallelToX1Axis()) return mLine->getX3Centroid()-mRad;
   else if(this->isParallelToX2Axis()) return mLine->getX3Centroid()-mRad; 
   else if(this->isParallelToX3Axis()) return mLine->getX3Minimum(); 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}	
/*==========================================================*/
double FeHalfDisc3D::getX3Maximum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX3Centroid()+mRad;
   else if(this->isParallelToX2Axis()) return mLine->getX3Centroid()+mRad; 
   else if(this->isParallelToX3Axis()) return mLine->getX3Maximum(); 
   else throw UbException(UB_EXARGS,"derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
bool FeHalfDisc3D::isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");
}
/*==========================================================*/
bool FeHalfDisc3D::isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");
}
/*=======================================================*/
bool FeHalfDisc3D::isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");
}

/*==========================================================*/
string FeHalfDisc3D::toString() 
{
	stringstream ss;
	ss<<"FeHalfDisc3D[";
	ss<<"line="<<this->mLine->toString();
   ss<<", r="<<this->mRad;
   ss<<"]";
   return(ss.str());
}
/*==========================================================*/
bool FeHalfDisc3D::isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");   
}
/*==========================================================*/
bool FeHalfDisc3D::isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");
}
/*==========================================================*/
GbLine3D* FeHalfDisc3D::createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2)
{
   throw UbException(UB_EXARGS,"sollte mal einer machen ... ");
}
/*==========================================================*/
vector<GbTriangle3D*> FeHalfDisc3D::getSurfaceTriangleSet()
{
   double x1ma,x1mb,x2m,x3m;
   if( this->isParallelToX1Axis() ) 
   {
      x1ma = this->getX1Minimum();
      x1mb = this->getX1Maximum();
      x2m  = this->getX2Centroid();
      x3m  = this->getX3Centroid();
   }
   else if( this->isParallelToX2Axis() ) 
   {
      x1ma = this->getX2Minimum();
      x1mb = this->getX2Maximum();
      x2m  = this->getX1Centroid();
      x3m  = this->getX3Centroid();
   }
   else if( this->isParallelToX3Axis() ) 
   {
      x1ma = this->getX3Minimum();
      x1mb = this->getX3Maximum();
      x2m  = this->getX2Centroid();
      x3m  = this->getX1Centroid();
   }
   else throw UbException(UB_EXARGS,"cylinder is not axis prallel");
   
   vector<GbTriangle3D*> triangles;    

   int segmentsCircle  = 14;
   double deltaPhi = UbMath::PI/(double)segmentsCircle;

   double phiX1a,phiX1b;
   double x1a,x2a,x3a,x1b,x2b,x3b,x1c,x2c,x3c,x1d,x2d,x3d;
   
   double dXCylinder =  fabs((x1mb-x1ma)); ///(double)0.5;
   int segmentsCylinder = (int)(fabs(x1mb-x1ma)/dXCylinder);
   for(int segCyl = 0; segCyl<segmentsCylinder; segCyl++)
   {
      x1a = x1d = x1ma+segCyl*dXCylinder;
      x1b = x1c = x1a+dXCylinder;
      
      for(phiX1a=UbMath::PI-deltaPhi; phiX1a>-deltaPhi; phiX1a-=deltaPhi)
      {
         phiX1b = phiX1a+deltaPhi;
         
         x2a =  x2m+mRad*std::sin(phiX1a);
         x3a =  x3m+mRad*std::cos(phiX1a);
         x2b =  x2m+mRad*std::sin(phiX1b);
         x3b =  x3m+mRad*std::cos(phiX1b);
         
         if( this->isParallelToX1Axis() ) 
         { 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1b,x2a,x3a),new FePoint3D(x1b,x2b,x3b),new FePoint3D(x1a,x2a,x3a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1a,x2b,x3b),new FePoint3D(x1a,x2a,x3a),new FePoint3D(x1b,x2b,x3b))); 
         }
         else if( this->isParallelToX2Axis() ) 
         { 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x2b,x1b,x3b),new FePoint3D(x2a,x1b,x3a),new FePoint3D(x2a,x1a,x3a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1a,x3a),new FePoint3D(x2b,x1a,x3b),new FePoint3D(x2b,x1b,x3b))); 
         }
         else if( this->isParallelToX3Axis() ) 
         { 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x3b,x2b,x1b),new FePoint3D(x3a,x2a,x1b),new FePoint3D(x3a,x2a,x1a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1a),new FePoint3D(x3b,x2b,x1a),new FePoint3D(x3b,x2b,x1b))); 
         }

      }
   }

   x2a = x2m;
   x3a = x3m;
   x2d = x2m;
   x3d = x3m+mRad;
   x3b = x3m;
   x3c = x3m-mRad;
   
   triangles.push_back(new GbTriangle3D(new FePoint3D(x1ma,x2a,x3a),new FePoint3D(x1mb,x2a,x3a),new FePoint3D(x1ma,x2a,x3d)));
   triangles.push_back(new GbTriangle3D(new FePoint3D(x1ma,x2a,x3d),new FePoint3D(x1mb,x2a,x3a),new FePoint3D(x1mb,x2a,x3d)));
   triangles.push_back(new GbTriangle3D(new FePoint3D(x1mb,x2a,x3b),new FePoint3D(x1ma,x2a,x3b),new FePoint3D(x1ma,x2a,x3c)));
   triangles.push_back(new GbTriangle3D(new FePoint3D(x1mb,x2a,x3b),new FePoint3D(x1ma,x2a,x3c),new FePoint3D(x1mb,x2a,x3c)));
  
   for(phiX1a=UbMath::PI-deltaPhi; phiX1a>-deltaPhi; phiX1a-=deltaPhi)
   {
      phiX1b = phiX1a+deltaPhi;

      x2a =  x2m;
      x3a =  x3m;
      x2b =  x2m;
      x3b =  x3m;
      x2c =  x2m+mRad*std::sin(phiX1b);
      x3c =  x3m+mRad*std::cos(phiX1b);
      x2d =  x2m+mRad*std::sin(phiX1a);
      x3d =  x3m+mRad*std::cos(phiX1a);

      if( this->isParallelToX1Axis() ) 
      { 
         triangles.push_back(new GbTriangle3D(new FePoint3D(x1ma,x2d,x3d),new FePoint3D(x1ma,x2c,x3c),new FePoint3D(x1ma,x2a,x3a)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x1mb,x2d,x3d),new FePoint3D(x1mb,x2a,x3a),new FePoint3D(x1mb,x2c,x3c)));
      }
      else if( this->isParallelToX2Axis() ) 
      { 
         triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1ma,x3a),new FePoint3D(x2b,x1ma,x3b),new FePoint3D(x2c,x1ma,x3c)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x2c,x1ma,x3c),new FePoint3D(x2d,x1ma,x3d),new FePoint3D(x2a,x1ma,x3a)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x2c,x1mb,x3c),new FePoint3D(x2b,x1mb,x3b),new FePoint3D(x2a,x1mb,x3a)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1mb,x3a),new FePoint3D(x2d,x1mb,x3d),new FePoint3D(x2c,x1mb,x3c)));
      }
      else if( this->isParallelToX3Axis() ) 
      { 
         triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1ma),new FePoint3D(x3b,x2b,x1ma),new FePoint3D(x3c,x2c,x1ma)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x3c,x2c,x1ma),new FePoint3D(x3d,x2d,x1ma),new FePoint3D(x3a,x2a,x1ma)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x3c,x2c,x1mb),new FePoint3D(x3b,x2b,x1mb),new FePoint3D(x3a,x2a,x1mb)));
         triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1mb),new FePoint3D(x3d,x2d,x1mb),new FePoint3D(x3c,x2c,x1mb)));
      }
   }

   return triangles;
}
/*==========================================================*/
void FeHalfDisc3D::objectChanged(UbObservable* changedObject)
{
   GbLine3D* line = dynamic_cast<GbLine3D*>(changedObject);
   if(!line || this->mLine!=line) return;

   this->notifyObserversObjectChanged();
}
/*==========================================================*/
void FeHalfDisc3D::objectWillBeDeleted(UbObservable* objectForDeletion)
{
   if(this->mLine)
   {
      UbObservable* observedObj = dynamic_cast<UbObservable*>(this->mLine);
      if(objectForDeletion == observedObj) { this->mLine = NULL; }
   }
}
/*=======================================================*/
void FeHalfDisc3D::scale(const double& sx1, const double& sx2, const double& sx3)
{  
   if( this->isParallelToX1Axis() )
   {
      if(!UbMath::equal(sx2,sx3)) throw UbException(UB_EXARGS,"|| to x1 -> different scaling sx2 and sx3 not possible");
      this->mRad*=sx2;
   }
   else if( this->isParallelToX2Axis() )
   {
      if(!UbMath::equal(sx1,sx3)) throw UbException(UB_EXARGS,"|| to x2 -> different scaling sx1 and sx3 not possible");
      this->mRad*=sx1;
   }
   else if( this->isParallelToX3Axis() )
   {
      if(!UbMath::equal(sx1,sx2)) throw UbException(UB_EXARGS,"|| to x3 -> different scaling sx1 and sx2 not possible");
      this->mRad*=sx1;
   }
   else throw UbException(UB_EXARGS,"unknown direction");

   this->mLine->scale(sx1,sx2,sx3);
   //notify observer wird automatisch aufgerufen
}
/*==========================================================*/
void FeHalfDisc3D::write(UbFileOutput* out) 
{                                      
   out->writeString(this->getCreator()->getTypeID());
   mLine->write(out);
   out->writeDouble(mRad);
   out->writeInteger(cylinderType);
}
/*==========================================================*/
void FeHalfDisc3D::read(UbFileInput* in) 
{  
   in->readString();                                    
   mLine = new GbLine3D;
   mLine->read(in);
   mRad         = in->readDouble();
   cylinderType = in->readInteger();
}
/*==========================================================*/
double FeHalfDisc3D::getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3)
{
   /*
   Distance D of the intersection between a Ray((ox1,ox2,ox3),(dx1,dx2,dx3)) and a Plane P: ax+by+cz+d=0 
   dc = a*dx1 + b*dx2 + c*dx3
   dw = a*ox1 + b*ox2 + c*ox3 + d
   D =   - dw / dc
   */
   double px1, px2, px3;
   double d = Ub::inf; // Distance to Min or Max Plane of the Cylinder  
                       // final distance should be less that d 
   
   if( this->isParallelToX1Axis() )
   {
      double minX1 = this->getX1Minimum();
      double maxX1 = this->getX1Maximum();

      if     (UbMath::equal(x1 ,minX1) && UbMath::negative(rx1))    return -1.0; 
      else if(UbMath::equal(x1 ,maxX1) && UbMath::positive(rx1))    return -1.0; 

      //falls die Linie nicht parallel zu den Seitenflächen ist
      if( x1< minX1  ||  x1 > maxX1 ) //nur für punkte links und rechts des cylinders
      {
         px1 = (x1 < minX1 ? minX1 : maxX1);
         //falls die Linie nicht parallel zu den Seitenflächen ist
         if( !UbMath::zero(rx1) )
         {
            // Plane a= 0, b= 1, c=0 d= -1*px2
            d   = -1.0*(x1 - px1) / rx1;
            px2 = x2 + d*rx2;
            px3 = x3 + d*rx3;

            if(UbMath::greater(mLine->getDistance(px1,px2,px3) , mRad))
            {
               if     (x1 < minX1 && rx1>0.0 ) d = Ub::inf;  //punkt liegt "links" vom cylinder und strahl hat evtl weiteren SP auf oberfläche 
               else if(x1 > maxX1 && rx1<0.0 ) d = Ub::inf;
               else return -1.0;
            }
            else return d;
         }
         else return -1.0;
      }
      else 
      {
         if     (UbMath::negative(rx1)) d = -1.0 * (x1 - minX1) / rx1;
         else if(UbMath::positive(rx1)) d = -1.0 * (x1 - maxX1) / rx1;
      }
   }
   else if( this->isParallelToX2Axis() )
   { 
      double minX2 = this->getX2Minimum();
      double maxX2 = this->getX2Maximum();

      if     (UbMath::equal(x2 ,minX2) && UbMath::negative(rx2))    return -1; 
      else if(UbMath::equal(x2 ,maxX2) && UbMath::positive(rx2))    return -1; 

      if( minX2 > x2  ||  x2 > maxX2 )
      {
         px2 = (x2 < minX2 ? minX2 : maxX2);
         //falls die Linie nicht parallel zu den Seitenflächen ist
         if( !UbMath::zero(rx2) )
         {
            // Plane a= 0, b= 1, c=0 d= -1*px2
            d   = -1*(x2 - px2) / rx2;
            px1 = x1 + d*rx1;
            px3 = x3 + d*rx3;
            
            if (UbMath::greater(mLine->getDistance(px1,px2,px3) , mRad))
            {
               if     (x2 < minX2 && rx2>0.0 ) d = Ub::inf;  //punkt liegt "links oberhalb" vom cylinder und strahl mit pos x1 hat evtl weiteren SP auf oberfläche 
               else if(x2 > maxX2 && rx2<0.0 ) d = Ub::inf;
               else return -1.0;
            }
            else return d;
         }
         else return -1.0;
      }
      else
      {
         if     (UbMath::negative(rx2)) d = -1.0 * (x2 - minX2) / rx2;
         else if(UbMath::positive(rx2)) d = -1.0 * (x2 - maxX2) / rx2;
      }
   }
   else if( this->isParallelToX3Axis() )
   {
      double minX3 = this->getX3Minimum();
      double maxX3 = this->getX3Maximum();

      if     (UbMath::equal(x3, minX3) && UbMath::negative(rx3)) return -1.0; 
      else if(UbMath::equal(x3, maxX3) && UbMath::positive(rx3)) return -1.0; 

      if(minX3 > x3  ||  x3 > maxX3 )
      {
         px3 = (x3 < minX3 ? minX3 : maxX3);
         //falls die Linie nicht parallel zu den Seitenflächen ist
         if (!UbMath::zero(rx3))
         {
            // Plane a= 0, b= 0, c=1 d= -1*px3
            d   = -1.0*(x3 - px3) / rx3;
            px2 = x2 + d*rx2;
            px1 = x1 + d*rx1;
            if( UbMath::greater(mLine->getDistance(px1,px2,px3) , mRad) )
            {
               if     (x3 < minX3 && rx3>0.0 ) d = Ub::inf;  
               else if(x3 > maxX3 && rx3<0.0 ) d = Ub::inf;
               else return -1.0;
            }
            else return d;
         }
         else return -1.0;
      }
      else 
      {
         if     (UbMath::negative(rx3)) d = -1.0 * (x3 - minX3) / rx3;
         else if(UbMath::positive(rx3)) d = -1.0 * (x3 - maxX3) / rx3;
      }
   }
   else throw UbException(UB_EXARGS,"funzt nur bei achsen parallelem cylinder");
   //////////////////////////////////////////////////////////////////////////
   //Q berechnen für Infinity Cylinder
   GbPoint3D* p1 = mLine->getPoint1();
   GbPoint3D* p2 = mLine->getPoint2();
   
   double axisX1 = p2->x1 - p1->x1;  /* Axis of the cylinder   */
   double axisX2 = p2->x2 - p1->x2;  /* mit p1 als base of cylinder */
   double axisX3 = p2->x3 - p1->x3;       

   //double dirlen = mLine->getLength(); 
   //double abs, t, s;

   double RCx1 = x1 - p1->x1;
   double RCx2 = x2 - p1->x2; 
   double RCx3 = x3 - p1->x3; 
   
   //n = ray x axis
   double nx1 = rx2*axisX3 - rx3*axisX2;
   double nx2 = rx3*axisX1 - rx1*axisX3;
   double nx3 = rx1*axisX2 - rx2*axisX1;
   double nLength = nx1*nx1 + nx2*nx2 + nx3*nx3;

   double abs;
   if( UbMath::zero( nLength ) )
   {  /* ray parallel to cyl  */
      //abs = RC dot axis
      double abs = RCx1*axisX1 + RCx2*axisX2 + RCx3*axisX3;
      double dx1 = RCx1 - abs*axisX1; 
      double dx2 = RCx2 - abs*axisX2;
      double dx3 = RCx3 - abs*axisX3;
      //abs   = sqrt(dx1*dx1 + dx2*dx2 + dx3*dx3);
      if( UbMath::greater( dx1*dx1 + dx2*dx2 + dx3*dx3 , mRad*mRad) ) 
         return -1.0;
   }

   //normalize "n"
   nLength = std::sqrt(nLength);
   double invnLength = 1.0/nLength;
   nx1*=invnLength;
   nx2*=invnLength;
   nx3*=invnLength;

   //shortest distance  = fabs( RC dot n )
   abs = fabs( RCx1*nx1 + RCx2*nx2 + RCx3*nx3 );     
   
   if( UbMath::lessEqual(abs, mRad) )
   {                    /* if ray hits cylinder */
      //Ox1 = RC x axis
      double Ox1 = RCx2*axisX3 - RCx3*axisX2;
      double Ox2 = RCx3*axisX1 - RCx1*axisX3;
      double Ox3 = RCx1*axisX2 - RCx2*axisX1;
      //t = - O dot n / nLength;
      double t = - (Ox1*nx1 + Ox2*nx2 + Ox3*nx3) / nLength;
      
      //O = n x axis;
      Ox1 = nx2*axisX3 - nx3*axisX2;
      Ox2 = nx3*axisX1 - nx1*axisX3;
      Ox3 = nx1*axisX2 - nx2*axisX1;

      //normalize O
      invnLength = 1.0/sqrt(Ox1*Ox1 + Ox2*Ox2 + Ox3*Ox3);
      Ox1*=invnLength;
      Ox2*=invnLength;
      Ox3*=invnLength;

      double s = fabs( sqrt(mRad*mRad - abs*abs) / (rx1*Ox1 + rx2*Ox2 + rx3*Ox3) );
      
      if( UbMath::greater(t-s,d) ) return -1.0;
      
      return  t - s;
   }
   
   return -1.0;
}
/*==========================================================*/
