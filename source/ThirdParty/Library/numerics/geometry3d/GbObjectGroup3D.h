//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBOBJECTGROUP3D_H
#define GBOBJECTGROUP3D_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF
#ifdef CAB_CTL
   #include <ctl.h>
#endif //CAB_CTL

#include <vector>
#include <cmath>

#include <numerics/geometry3d/GbObject3D.h>
#include <basics/utilities/UbObserver.h>
#include <numerics/geometry3d/GbPoint3D.h>

#include <PointerDefinitions.h>

class GbLine3D;
class GbTriangle3D;
class GbObject3DCreator;

class GbObjectGroup3D : public GbObject3D, public UbObserver
{                                              
public:
   enum TRIANGULATIONMODE { CUBOIDPROJECTION ,RAYPROJECTION };
   
   //////////////////////////////////////////////////////////////////////////
   // Konstruktoren
   GbObjectGroup3D(); 
   GbObjectGroup3D(GbObjectGroup3D *group){}; 
   ~GbObjectGroup3D();

   GbObjectGroup3D* clone() {return new GbObjectGroup3D(this);}
   void finalize();

   void addGbObject(GbObject3D* object)
   {
      this->geoobjects.push_back(object);
   }

   double getRadius() const	{	return this->radius;	}

   double getX1Centroid()  { return midPoint->getX1Coordinate();}
   double getX1Minimum()   { return midPoint->getX1Coordinate()-radius;}
   double getX1Maximum()   { return midPoint->getX1Coordinate()+radius;}
   double getX2Centroid()  { return midPoint->getX2Coordinate();}
   double getX2Minimum()   { return midPoint->getX2Coordinate()-radius;}
   double getX2Maximum()   { return midPoint->getX2Coordinate()+radius;}
   double getX3Centroid()  { return midPoint->getX3Coordinate();}
   double getX3Minimum()   { return midPoint->getX3Coordinate()-radius;}
   double getX3Maximum()   { return midPoint->getX3Coordinate()+radius;}

   void setCenterX1Coordinate(const double& value);
   void setCenterX2Coordinate(const double& value);
   void setCenterX3Coordinate(const double& value);
   void setCenterCoordinates(const double& x1, const double& x2, const double& x3);
   void setRadius(const double& radius);

   GbLine3D* createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2);
   double getDistance(GbPoint3D* p); 
   double getDistance(const double& x1p, const double& x2p, const double& x3p);

   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3);
   bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary);

   bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   double getCellVolumeInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   
   std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   void addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles);
                                
   bool hasRaytracing() { return true; }
   /*|r| must be 1! einheitsvector!!*/
   double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3);

   bool hasIntersectionWithDirectedLine(GbPoint3D origin, GbPoint3D direction);

	std::string toString();

   ObObjectCreator* getCreator(){ return NULL; };
   void write(UbFileOutput* out);
   void read(UbFileInput* in);       

   void translate(const double& x1, const double& x2, const double& x3) 
   {
      this->midPoint->translate(x1, x2, x3); 
      this->notifyObserversObjectChanged();
   }
   void rotate(const double& rx1, const double& rx2, const double& rx3) {/* rotation makes no sense*/ }
   void scale(const double& sx1, const double& sx2, const double& sx3) { this->radius *= sx1; }

   TRIANGULATIONMODE getTriangulationMode() {return triangulationMode;}
   void setTriangulationMode(TRIANGULATIONMODE mode) { this->triangulationMode = mode; }
   
   //virtuelle Methoden von UbObserver
   void objectChanged(UbObservable* changedObject)
   {
      this->notifyObserversObjectChanged();
      //std::cout<<"GbSphere:objectChanged() - toDo-);";
   }
   void objectWillBeDeleted(UbObservable* objectForDeletion)
   {
	   std::cout<<"throw UbException(-GbObjectGroup3D::finalize() - toDo-);";
   }

   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere, weil man eine

   std::list< GbObject3D* > getGbObject3DList() { return this->geoobjects; }

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & midPoint;
      ar & radius;
      ar & triangulationMode;
      ar & geoobjects;
   }
#endif //CAB_RCF
#ifdef CAB_CTL
   ctl::oStream &write(ctl::oStream &os) const
   { 
      midPoint->write(os);
      return os<<radius; 
   }
   ctl::iStream &read(ctl::iStream &is) 
   { 
      midPoint->read(is);
      return is>>radius;
   }
#endif //CAB_CTL

private:
   GbPoint3D* midPoint;
   double radius;  // Radius des Kreises
   TRIANGULATIONMODE triangulationMode;

   std::list< GbObject3D* > geoobjects;
};

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   #if defined(RCF_USE_SF_SERIALIZATION) && (CAB_RCF <= 903) 
      SF_SERIALIZE_ENUM(GbObjectGroup3D::TRIANGULATIONMODE) //bei klassen ausserhalb der klasse;-)
   #endif
   UB_AUTO_RUN_NAMED(   SF::registerType<GbObjectGroup3D>("GbObjectGroup3D")        , SF_GbObjectGroup3D     );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbObjectGroup3D >()), SF_GbObjectGroup3D_BD1 );
#endif //RCF_USE_SF_SERIALIZATION

#endif //GbObjectGroup3D_H
