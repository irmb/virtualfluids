#include <numerics/geometry3d/fem/FeRing3D.h>    
#include <numerics/geometry3d/GbSystem3D.h>
#include <numerics/geometry3d/fem/FePoint3D.h>                    
#include <numerics/geometry3d/GbLine3D.h>
#include <numerics/geometry3d/GbTriangle3D.h>

using namespace std;

/*=======================================================*/
ObObjectCreator* FeRing3D::getCreator()
{
   return NULL; //FeRing3DCreator::getInstance(); 
}
// Konstruktor
/*==========================================================*/
FeRing3D::FeRing3D()
{
   GbPoint3D* p1 = new GbPoint3D();
   GbPoint3D* p2 = new GbPoint3D();
   mLine = new GbLine3D(p1,p2);
   this->mLine->addObserver(this);
   inRadius = 0.0;
   outRadius = 0.0;
   ringType = FeRing3D::NOTPARALLELTOAXIS;
}                                                   
/*=======================================================*/
FeRing3D::FeRing3D(FeRing3D* ring)
{
   inRadius     = ring->getInRadius();
   outRadius    = ring->getOutRadius();
   ringType     = ring->ringType;
   mLine        = ring->getLine()->clone();

   this->mLine->addObserver(this);
}                                                   
/*==========================================================*/
FeRing3D::FeRing3D(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b, const double& inradius, const double& outradius)
{
   mLine = new GbLine3D;
   mLine->setPoints( new GbPoint3D(min(x1a,x1b), min(x2a,x2b), min(x3a,x3b))
	                 ,new GbPoint3D(max(x1a,x1b), max(x2a,x2b), max(x3a,x3b)));
   this->mLine->addObserver(this);
   this->inRadius = inradius;
   this->outRadius = outradius;

   this->initRingType();
   if((this->ringType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException("FeRing3D::FeRing3D - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}                                                            
/*==========================================================*/
FeRing3D::FeRing3D(GbPoint3D* p1, GbPoint3D* p2, const double& inradius, const double& outradius)
{
   this->inRadius = inradius;
   this->outRadius = outradius;

   mLine = new GbLine3D(p1,p2);
   this->mLine->addObserver(this);
   this->initRingType();
   if((this->ringType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException("FeRing3D::FeRing3D - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
FeRing3D::FeRing3D(GbLine3D* line, const double& inradius, const double& outradius)
{
   this->inRadius = inradius;
   this->outRadius = outradius;

   this->mLine = line;
   this->mLine->addObserver(this);
   this->initRingType();
   if((this->ringType & NOTPARALLELTOAXIS)==NOTPARALLELTOAXIS) 
      throw UbException("FeRing3D::FeRing3D - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/		
// Destruktor
FeRing3D::~FeRing3D()
{
   if(mLine) this->mLine->removeObserver(this);
   mLine = NULL;
}
/*=======================================================*/
void FeRing3D::initRingType()
{
   double x1a = mLine->getPoint1()->x1;    double x1b = mLine->getPoint2()->x1;
   double x2a = mLine->getPoint1()->x2;    double x2b = mLine->getPoint2()->x2;
   double x3a = mLine->getPoint1()->x3;    double x3b = mLine->getPoint2()->x3;
   
   if     (x1a!=x1b && x2a==x2b && x3a==x3b)  this->ringType = X1PARALLEL; 
   else if(x2a!=x2b && x1a==x1b && x3a==x3b)  this->ringType = X2PARALLEL; 
   else if(x3a!=x3b && x1a==x1b && x2a==x2b)  this->ringType = X3PARALLEL; 
   else                                       this->ringType = NOTPARALLELTOAXIS;
}
/*=======================================================*/
void FeRing3D::finalize() 
{ 
   if(this->mLine) 
   {
      mLine->finalize();
      delete mLine; 
      mLine=NULL;
   } 
}
/*=======================================================*/
double FeRing3D::getHeight()
{
   if(mLine) return mLine->getLength(); return 0.0; 
}
/*=======================================================*/
GbPoint3D* FeRing3D::getPoint1()
{
   if(this->mLine) return this->mLine->getPoint1();
   return NULL;
}
/*=======================================================*/
GbPoint3D* FeRing3D::getPoint2()
{
   if(this->mLine) return this->mLine->getPoint2();
   return NULL;
}
/*=======================================================*/
void FeRing3D::setInRadius(const double& radius) 
{ 
   this->inRadius = std::fabs(radius); 
   this->notifyObserversObjectChanged();
}
/*=======================================================*/
void FeRing3D::setOutRadius(const double& radius) 
{ 
   this->outRadius = std::fabs(radius); 
   this->notifyObserversObjectChanged();
}
/*=======================================================*/
void FeRing3D::setLine(GbLine3D* line) 
{
   if(this->mLine) this->mLine->removeObserver(this);
   this->mLine = line;  
   this->mLine->addObserver(this);
   this->initRingType();

   this->notifyObserversObjectChanged();
}
/*=======================================================*/
void FeRing3D::setPoint1(const double& x1, const double& x2, const double& x3)
{ 
   if(!mLine->getPoint1()) throw UbException("FeRing3D::setPoint1() - line has no point1");
   mLine->getPoint1()->setCoordinates(x1,x2,x3);
   this->initRingType();

   //this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*=======================================================*/
void FeRing3D::setPoint2(const double& x1, const double& x2, const double& x3)
{ 
   if(!mLine->getPoint2()) throw UbException("FeRing3D::setPoint1() - line has no point2");
   mLine->getPoint2()->setCoordinates(x1,x2,x3);
   this->initRingType();

   //this->notifyObserversObjectChanged(); //wird automatisch aufgerufen, da der point (this) benachrichtigt...
}
/*==========================================================*/		
double FeRing3D::getX1Centroid()  
{
   return mLine->getX1Centroid();
}
/*==========================================================*/
double FeRing3D::getX1Minimum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX1Minimum(); 
   else if(this->isParallelToX2Axis()) return mLine->getX1Centroid()-outRadius; 
   else if(this->isParallelToX3Axis()) return mLine->getX1Centroid()-outRadius; 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Minimum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeRing3D::getX1Maximum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX1Maximum(); 
   else if(this->isParallelToX2Axis()) return mLine->getX1Centroid()+outRadius; 
   else if(this->isParallelToX3Axis()) return mLine->getX1Centroid()+outRadius; 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Maximum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeRing3D::getX2Centroid()
{
   return mLine->getX2Centroid();
}
/*==========================================================*/
double FeRing3D::getX2Minimum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX2Centroid()-outRadius;
   else if(this->isParallelToX2Axis()) return mLine->getX2Minimum();
   else if(this->isParallelToX3Axis()) return mLine->getX2Centroid()-outRadius; 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Minimum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}	
/*==========================================================*/
double FeRing3D::getX2Maximum()   
{
   if     (this->isParallelToX1Axis())  return mLine->getX2Centroid()+outRadius;
   else if(this->isParallelToX2Axis())  return mLine->getX2Maximum();
   else if(this->isParallelToX3Axis())  return mLine->getX2Centroid()+outRadius; 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Maximum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
double FeRing3D::getX3Centroid()
{
   return mLine->getX3Centroid();
}
/*==========================================================*/
double FeRing3D::getX3Minimum()   
{	
   if     (this->isParallelToX1Axis()) return mLine->getX3Centroid()-outRadius;
   else if(this->isParallelToX2Axis()) return mLine->getX3Centroid()-outRadius; 
   else if(this->isParallelToX3Axis()) return mLine->getX3Minimum(); 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Minimum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}	
/*==========================================================*/
double FeRing3D::getX3Maximum()   
{
   if     (this->isParallelToX1Axis()) return mLine->getX3Centroid()+outRadius;
   else if(this->isParallelToX2Axis()) return mLine->getX3Centroid()+outRadius; 
   else if(this->isParallelToX3Axis()) return mLine->getX3Maximum(); 
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getX3Maximum - derzeit nur zu Achsen orthogonale Cylinder erlaubt... isPointInObject3D funzt sonst ned");
}
/*==========================================================*/
bool FeRing3D::isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");

}
/*==========================================================*/
bool FeRing3D::isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}
/*=======================================================*/
bool FeRing3D::isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}

/*==========================================================*/
string FeRing3D::toString() 
{
	stringstream ss;
	ss<<"FeRing3D[";
	ss<<"line="<<this->mLine->toString();
   ss<<", inRadius="<<this->inRadius;
   ss<<", outRadius="<<this->outRadius;
   ss<<"]";
   return(ss.str());
}
/*==========================================================*/
bool FeRing3D::isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}
/*==========================================================*/
bool FeRing3D::isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}
/*==========================================================*/
GbLine3D* FeRing3D::createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}
/*==========================================================*/
vector<GbTriangle3D*> FeRing3D::getSurfaceTriangleSet()
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
   else throw UbException(__FILE__, __LINE__, "FeRing3D::getSurfaceTriangleSet() - ring not axis prallel");
   
   vector<GbTriangle3D*> triangles;    

   int segmentsCircle  = 14;
   double deltaPhi = UbMath::PI/(double)segmentsCircle;

   double phiX1a,phiX1b;
   double x1a,x2a,x3a,x1b,x2b,x3b,x1c,x2c,x3c,x1d,x2d,x3d;
   double x2aa,x3aa,x2bb,x3bb;
   
   double dXCylinder =  fabs((x1mb-x1ma)); // /(double)segmentsCircle;
   int segmentsCylinder = (int)(fabs(x1mb-x1ma)/dXCylinder);
   for(int segCyl = 0; segCyl<segmentsCylinder; segCyl++)
   {
      x1a = x1d = x1ma+segCyl*dXCylinder;
      x1b = x1c = x1a+dXCylinder;
      
      for(phiX1a=2.0*UbMath::PI; phiX1a>0; phiX1a-=deltaPhi)
      {
         phiX1b = phiX1a+deltaPhi;
         
         x2a =  x2m+this->outRadius*std::sin(phiX1a);
         x3a =  x3m+this->outRadius*std::cos(phiX1a);
         x2b =  x2m+this->outRadius*std::sin(phiX1b);
         x3b =  x3m+this->outRadius*std::cos(phiX1b);

         x2aa =  x2m+this->inRadius*std::sin(phiX1a);
         x3aa =  x3m+this->inRadius*std::cos(phiX1a);
         x2bb =  x2m+this->inRadius*std::sin(phiX1b);
         x3bb =  x3m+this->inRadius*std::cos(phiX1b);
         
         if( this->isParallelToX1Axis() ) 
         { 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1b,x2a,x3a),new FePoint3D(x1b,x2b,x3b),new FePoint3D(x1a,x2a,x3a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1a,x2b,x3b),new FePoint3D(x1a,x2a,x3a),new FePoint3D(x1b,x2b,x3b))); 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1b,x2bb,x3bb),new FePoint3D(x1b,x2aa,x3aa),new FePoint3D(x1a,x2aa,x3aa)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1a,x2aa,x3aa),new FePoint3D(x1a,x2bb,x3bb),new FePoint3D(x1b,x2bb,x3bb))); 
         }
         else if( this->isParallelToX2Axis() ) 
         { 
            throw UbException(__FILE__,__LINE__," sollte man mal machen");
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2b,x1b,x3b),new FePoint3D(x2a,x1b,x3a),new FePoint3D(x2a,x1a,x3a)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1a,x3a),new FePoint3D(x2b,x1a,x3b),new FePoint3D(x2b,x1b,x3b))); 
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2bb,x1b,x3bb),new FePoint3D(x2aa,x1b,x3aa),new FePoint3D(x2aa,x1a,x3aa)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2aa,x1a,x3aa),new FePoint3D(x2bb,x1a,x3bb),new FePoint3D(x2bb,x1b,x3bb))); 
         }
         else if( this->isParallelToX3Axis() ) 
         { 
            throw UbException(__FILE__,__LINE__," sollte man mal machen");
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3b,x2b,x1b),new FePoint3D(x3a,x2a,x1b),new FePoint3D(x3a,x2a,x1a)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1a),new FePoint3D(x3b,x2b,x1a),new FePoint3D(x3b,x2b,x1b))); 
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3bb,x2bb,x1b),new FePoint3D(x3aa,x2aa,x1b),new FePoint3D(x3aa,x2aa,x1a)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3aa,x2aa,x1a),new FePoint3D(x3bb,x2bb,x1a),new FePoint3D(x3bb,x2bb,x1b))); 
         }

      }
   }
   
   //int segmentsSide = (int)(this->outRadius/dXCylinder);

   double radius0, radius1;
   //for(int segCyl = 10; segCyl<segmentsSide; segCyl++)
   //{
      radius0 = inRadius;//segCyl*dXCylinder;
      radius1 = outRadius;//radius0+dXCylinder;
      //if(segCyl==segmentsSide-1) radius1=outRadius;

      for(phiX1a=2.0*UbMath::PI; phiX1a>0; phiX1a-=deltaPhi)
      {
         phiX1b = phiX1a+deltaPhi;

         x2a =  x2m+radius0*std::sin(phiX1a);
         x3a =  x3m+radius0*std::cos(phiX1a);
         x2b =  x2m+radius0*std::sin(phiX1b);
         x3b =  x3m+radius0*std::cos(phiX1b);
         x2c =  x2m+radius1*std::sin(phiX1b);
         x3c =  x3m+radius1*std::cos(phiX1b);
         x2d =  x2m+radius1*std::sin(phiX1a);
         x3d =  x3m+radius1*std::cos(phiX1a);

         if( this->isParallelToX1Axis() ) 
         { 
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1ma,x2b,x3b),new FePoint3D(x1ma,x2a,x3a),new FePoint3D(x1ma,x2c,x3c)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1ma,x2d,x3d),new FePoint3D(x1ma,x2c,x3c),new FePoint3D(x1ma,x2a,x3a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1mb,x2b,x3b),new FePoint3D(x1mb,x2c,x3c),new FePoint3D(x1mb,x2a,x3a)));
            triangles.push_back(new GbTriangle3D(new FePoint3D(x1mb,x2d,x3d),new FePoint3D(x1mb,x2a,x3a),new FePoint3D(x1mb,x2c,x3c)));
         }                                                                   
         else if( this->isParallelToX2Axis() ) 
         { 
            throw UbException(__FILE__,__LINE__," sollte man mal machen");
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1ma,x3a),new FePoint3D(x2b,x1ma,x3b),new FePoint3D(x2c,x1ma,x3c)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2c,x1ma,x3c),new FePoint3D(x2d,x1ma,x3d),new FePoint3D(x2a,x1ma,x3a)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2c,x1mb,x3c),new FePoint3D(x2b,x1mb,x3b),new FePoint3D(x2a,x1mb,x3a)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x2a,x1mb,x3a),new FePoint3D(x2d,x1mb,x3d),new FePoint3D(x2c,x1mb,x3c)));
         }
         else if( this->isParallelToX3Axis() ) 
         { 
            throw UbException(__FILE__,__LINE__," sollte man mal machen");
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1ma),new FePoint3D(x3b,x2b,x1ma),new FePoint3D(x3c,x2c,x1ma)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3c,x2c,x1ma),new FePoint3D(x3d,x2d,x1ma),new FePoint3D(x3a,x2a,x1ma)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3c,x2c,x1mb),new FePoint3D(x3b,x2b,x1mb),new FePoint3D(x3a,x2a,x1mb)));
            //triangles.push_back(new GbTriangle3D(new FePoint3D(x3a,x2a,x1mb),new FePoint3D(x3d,x2d,x1mb),new FePoint3D(x3c,x2c,x1mb)));
         }
      }
  // }

   return triangles;
}
/*==========================================================*/
void FeRing3D::objectChanged(UbObservable* changedObject)
{
   GbLine3D* line = dynamic_cast<GbLine3D*>(changedObject);
   if(!line || this->mLine!=line) return;

   this->notifyObserversObjectChanged();
}
/*==========================================================*/
void FeRing3D::objectWillBeDeleted(UbObservable* objectForDeletion)
{
   if(this->mLine)
   {
      UbObservable* observedObj = dynamic_cast<UbObservable*>(this->mLine);
      if(objectForDeletion == observedObj) { this->mLine = NULL; }
   }
}
/*=======================================================*/
void FeRing3D::scale(const double& sx1, const double& sx2, const double& sx3)
{  
   if( this->isParallelToX1Axis() )
   {
      if(!UbMath::equal(sx2,sx3)) throw UbException("FeRing3D::scale - || to x1 -> different scaling sx2 and sx3 not possible");
      this->inRadius*=sx2;
      this->outRadius*=sx2;
   }
   else if( this->isParallelToX2Axis() )
   {
      if(!UbMath::equal(sx1,sx3)) throw UbException("FeRing3D::scale - || to x2 -> different scaling sx1 and sx3 not possible");
      this->inRadius*=sx1;
      this->outRadius*=sx1;
   }
   else if( this->isParallelToX3Axis() )
   {
      if(!UbMath::equal(sx1,sx2)) throw UbException("FeRing3D::scale - || to x3 -> different scaling sx1 and sx2 not possible");
      this->inRadius*=sx1;
      this->outRadius*=sx1;
   }
   else throw UbException("FeRing3D::scale - unknown direction");

   this->mLine->scale(sx1,sx2,sx3);
   //notify observer wird automatisch aufgerufen
}
/*==========================================================*/
void FeRing3D::write(UbFileOutput* out) 
{                                      
   out->writeString(this->getCreator()->getTypeID());
   mLine->write(out);
   out->writeDouble(inRadius);
   out->writeDouble(outRadius);
   out->writeInteger(ringType);
}
/*==========================================================*/
void FeRing3D::read(UbFileInput* in) 
{  
   in->readString();                                    
   mLine = new GbLine3D;
   mLine->read(in);
   inRadius  = in->readDouble();
   outRadius = in->readDouble();
   ringType  = in->readInteger();
}
/*==========================================================*/
double FeRing3D::getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3)
{
   throw UbException(__FILE__,__LINE__,"FeRing3D function not implemented");
}
/*==========================================================*/

