//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBPOINT3D_H
#define GBPOINT3D_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <string>
#include <sstream>
#include <cmath>

#include <numerics/geometry3d/GbObject3D.h>

#include <PointerDefinitions.h>

class GbTriangle3D;
class GbObject3DCreator;

#ifdef CAB_CTL
   #include <ctl.h>
#endif

class GbPoint3D : public GbObject3D
{
public:
   GbPoint3D();
   GbPoint3D(const double& x1, const double& x2, const double& x3);
   GbPoint3D(GbPoint3D *point);                
   ~GbPoint3D() {}

   GbPoint3D* clone() {return new GbPoint3D(this);}
   void finalize() {}

   void setCoordinates(const double& x1, const double& x2, const double& x3)
   {
       this->x1=x1;
       this->x2=x2;
       this->x3=x3;
       this->notifyObserversObjectChanged();
   }
   void setX1(const double& x1) { this->x1=x1; this->notifyObserversObjectChanged(); }
   void setX2(const double& x2) { this->x2=x2; this->notifyObserversObjectChanged(); }
   void setX3(const double& x3) { this->x3=x3; this->notifyObserversObjectChanged(); }

   double getX1Coordinate() const  { return this->x1; }
   double getX2Coordinate() const  { return this->x2; }
   double getX3Coordinate() const  { return this->x3; }

   void transform(const double matrix[4][4]);
 
   double getX1Centroid()  { return this->x1; }
   double getX1Minimum()   { return this->x1; }
   double getX1Maximum()   { return this->x1; }
   double getX2Centroid()  { return this->x2; }
   double getX2Minimum()   { return this->x2; }
   double getX2Maximum()   { return this->x2; }
   double getX3Centroid()  { return this->x3; }
   double getX3Minimum()   { return this->x3; }
   double getX3Maximum()   { return this->x3; }        
 
   void translate(const double& x1, const double& x2, const double& x3);
   void rotate(const double& rx1, const double& rx2, const double& rx3);
   void scale(const double& sx1, const double& sx2, const double& sx3);

   double getDistance(GbPoint3D *p);
   bool equals(const GbPoint3D* point) const;
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary);
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3);
   bool isCellInsideGbObject3D(const double& x11,const double& x21,const double& x31,const double& x12,const double& x22,const double& x23) { return false; }

   std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   GbLine3D* createClippedLine3D(GbPoint3D &point1, GbPoint3D &point2);
   virtual std::string toString();
   ObObjectCreator* getCreator();
   void write(UbFileOutput* out);
   void read(UbFileInput* in);

   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren
                                          //, welche sonst hier "ueberdeckt" waere,da es dieselbe methode mit anderen args gibt!
#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & x1; 
      ar & x2; 
      ar & x3;
   }
#endif //CAB_RCF

#ifdef CAB_CTL
   ctl::oStream &write(ctl::oStream &os) const
   { 
      return os<<x1<<x2<<x3; 
   }
   ctl::iStream &read(ctl::iStream &is) 
   { 
      return is>>x1>>x2>>x3;
   }
#endif

   //member
   double x1;
   double x2;
   double x3;      
};


#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   UB_AUTO_RUN_NAMED(   SF::registerType<GbPoint3D>("GbPoint3D")              , SF_GbPoint3D      );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbPoint3D >()), SF_GbPoint3D_BD1 );
#endif //RCF_USE_SF_SERIALIZATION

#endif
