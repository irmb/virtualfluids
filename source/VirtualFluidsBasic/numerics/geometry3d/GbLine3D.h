//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBLINE3D_H
#define GBLINE3D_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <sstream>
#include <cmath>
          
#include <basics/utilities/UbObserver.h>

#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/GbPoint3D.h>

class GbCuboid3D;
class GbObject3DCreator;

#include <PointerDefinitions.h>


/*=========================================================================*/
/* GbLine3D                                                                */
/*                                                                         */
/**
 * This Class provides basic 3D line objects.
 * The describing points are observed by 3D line objects.
 * <BR><BR><HR>
*/

class GbLine3D	: public GbObject3D , public UbObserver
{
public:
   GbLine3D();
	GbLine3D(GbPoint3D* point1, GbPoint3D* point2);
	GbLine3D(GbLine3D* line);
   ~GbLine3D(); 

   GbLine3D* clone() { return new GbLine3D(this); }
   void finalize();

   void setPoint1(GbPoint3D* point1);
   void setPoint2(GbPoint3D* point2);
   void setPoints(GbPoint3D* point1, GbPoint3D* point2);

   void deletePoint1() { if(this->p1) {this->p1->removeObserver(this); delete this->p1; this->p1=NULL;} }
   void deletePoint2() { if(this->p2) {this->p2->removeObserver(this); delete this->p2; this->p2=NULL;} }
   void deletePoints() { this->deletePoint1(); this->deletePoint2(); }

   GbPoint3D* getPoint1() { return this->p1; }
   GbPoint3D* getPoint2() { return this->p2; }    
   
   double getLength()     { return(this->length); }
	
   double getX1Centroid() { return((this->p1->x1+this->p2->x1)*0.5);}
   double getX2Centroid() { return((this->p1->x2+this->p2->x2)*0.5); };
   double getX3Centroid() { return((this->p1->x3+this->p2->x3)*0.5); }
   
   double getX1Minimum()  { return(this->p1->x1 < this->p2->x1 ? this->p1->x1 : this->p2->x1); }
   double getX2Minimum()  { return(this->p1->x2 < this->p2->x2 ? this->p1->x2 : this->p2->x2); }
   double getX3Minimum()  { return(this->p1->x3 < this->p2->x3 ? this->p1->x3 : this->p2->x3); }
   
   double getX1Maximum()  { return(this->p1->x1 > this->p2->x1 ? this->p1->x1 : this->p2->x1); }
   double getX2Maximum()  { return(this->p1->x2 > this->p2->x2 ? this->p1->x2 : this->p2->x2); }
   double getX3Maximum()  { return(this->p1->x3 > this->p2->x3 ? this->p1->x3 : this->p2->x3); }
	                                               
   void scale(const double& sx1, const double& sx2, const double& sx3);
   void translate(const double& tx1, const double& tx2, const double& tx3);

   GbPoint3D* calculateIntersectionPoint3D(GbLine3D* line);
   GbLine3D*  createClippedLine3D(GbCuboid3D* cuboid);
   GbLine3D*  createClippedLine3D(GbPoint3D* pA, GbPoint3D* pE);
   
   double     getDistance(const GbPoint3D& point);
   double     getDistance(const double& x1,const double& x2,const double& x3);

   std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
   {
      throw UbException(UB_EXARGS,"not implemented");
   }
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary)
   {
      throw UbException(UB_EXARGS,"not implemented");
   }
   bool isCellInsideGbObject3D(const double& x11,const double& x21,const double& x31,const double& x12,const double& x22,const double& x32) { return false; }

   GbLine3D* createClippedLine3D (GbPoint3D& point1, GbPoint3D& point2)
   {
      throw UbException(UB_EXARGS,"not implemented");
   }

   //virtuelle Methoden von UbObserver
   void objectChanged(UbObservable* changedObject);
   void objectWillBeDeleted(UbObservable* objectForDeletion);

   std::string toString();
   ObObjectCreator* getCreator();
   void write(UbFileOutput* out);
   void read(UbFileInput* in);

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & p1;
      ar & p2;
      ar & length;
      if( ArchiveTools::isReading(ar) ) this->calculateValues();
   }
#endif //CAB_RCF

   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere
protected:
   GbPoint3D* p1;
	GbPoint3D* p2;
	double     length;

private:
   void calculateValues();
};

#ifdef RCF_USE_SF_SERIALIZATION
    UB_AUTO_RUN_NAMED(   SF::registerType<GbLine3D>("GbLine3D"), SF_GbLine3D  );
    UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbLine3D >()), SF_GbLine3D_BD1 );
    UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< UbObserver, GbLine3D>() ), SF_GbLine3D_BD2 );
#endif //RCF_USE_SF_SERIALIZATION

#endif
