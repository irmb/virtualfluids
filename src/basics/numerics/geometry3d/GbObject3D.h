//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBOBJECT3D_H
#define GBOBJECT3D_H

#include <string>
#include <vector>


#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <basics/utilities/UbSystem.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>
#include <basics/utilities/UbFileOutput.h>
#include <basics/utilities/UbObservable.h>
#include <basics/utilities/UbTuple.h>
#include <basics/objects/ObObject.h>

class GbPoint3D;
class GbLine3D;
class GbTriangle3D;
class GbObject3DCreator;

#ifdef CAB_CTL
#include <ctl.h>
#endif

#include "basics_export.h"
#include <PointerDefinitions.h>


/*=========================================================================*/
/* GbObject3D                                                              */
/*                                                                         */
/**
 * This Interface provides basic 3D geometry objects methods.
 * <BR><BR><HR>
 * @author <A HREF="mailto:geller@cab.bau.tu-bs.de">S. Geller</A>
 * @author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
 * @version 1.0 - 02.02.05
*/
class BASICS_EXPORT GbObject3D : public ObObject
{
public:
#ifdef CAB_CTL
   virtual ctl::oStream &write(ctl::oStream &os) const
   {
      return os;
   }
   virtual ctl::iStream &read(ctl::iStream &is)
   {
      return is;
   }
#endif

   virtual ~GbObject3D(){}

   //ueberschriebene methode von ObObject
   virtual std::string getTypeID();

   //abstract Methods
   virtual void finalize() =0 ; //detroys also all dynamic objects (e.g. GbPoints in GbLine)
   virtual ObObjectCreator* getCreator()=0;

   /**
    * Returns the centroid x1 coordinate of this 3D object.
    * @return the centroid x1 coordinate of this 3D object
    */
   virtual double getX1Centroid()=0;
   /**
    * Returns the minimum x1 coordinate of this 3D object.
    * @return the minimum x1 coordinate of this 3D object
    */
   virtual double getX1Minimum()=0;
   /**
    * Returns the maximum x1 coordinate of this 3D object.
    * @return the maximum x1 coordinate of this 3D object
    */
   virtual double getX1Maximum()=0;
   /**
    * Returns the centroid x2 coordinate of this 3D object.
    * @return the centroid x2 coordinate of this 3D object
    */
   virtual double getX2Centroid()=0;
   /**
    * Returns the minimum x2 coordinate of this 3D object.
    * @return the minimum x2 coordinate of this 3D object
    */
   virtual double getX2Minimum()=0;
   /**
    * Returns the maximum x2 coordinate of this 3D object.
    * @return the maximum x2 coordinate of this 3D object
    */
   virtual double getX2Maximum()=0;

	virtual double getX3Centroid()=0;
   /**
    * Returns the minimum x2 coordinate of this 3D object.
    * @return the minimum x2 coordinate of this 3D object
    */
   virtual double getX3Minimum()=0;
   /**
    * Returns the maximum x2 coordinate of this 3D object.
    * @return the maximum x2 coordinate of this 3D object
    */
   virtual double getX3Maximum()=0;

   /*=======================================================*/
   double getLengthX1() { return (getX1Maximum()-getX1Minimum()); }
   double getLengthX2() { return (getX2Maximum()-getX2Minimum()); }
   double getLengthX3() { return (getX3Maximum()-getX3Minimum()); }

   virtual void setCenterX1Coordinate(const double& value) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }
   virtual void setCenterX2Coordinate(const double& value) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }
   virtual void setCenterX3Coordinate(const double& value) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }
   virtual void setCenterCoordinates(const double& x1, const double& x2, const double& x3) { throw UbException(UB_EXARGS, "not implemented for " + (std::string)typeid(*this).name()); }
   virtual void setCenterCoordinates(const UbTupleDouble3& position) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }

   //Rotates the Point in relation to the origen.
   //Parameters must be radian measure.
   virtual void rotate(const double& rx1, const double& rx2, const double& rx3) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }
   virtual void translate(const double& x1, const double& x2, const double& x3) { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }
   virtual void scale(const double& sx1, const double& sx2, const double& sx3)  { throw UbException(UB_EXARGS,"not implemented for "+(std::string)typeid(*this).name() ); }

   virtual void write(UbFileOutput* out)=0;
   virtual void read(UbFileInput* in)=0;

   virtual bool isPointInGbObject3D(GbPoint3D* p);
   virtual bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary)=0;
   virtual bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3)=0;

   virtual bool isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   virtual bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   virtual bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   virtual double getCellVolumeInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b){ return -1.0;};

   virtual bool isInsideCell(const double& minX1,const double& minX2,const double& minX3,const double& maxX1,const double& maxX2,const double& maxX3);

   virtual GbLine3D* createClippedLine3D (GbPoint3D &point1, GbPoint3D &point2)=0;
   virtual std::vector<GbTriangle3D*> getSurfaceTriangleSet()=0;

   virtual void addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles) { throw UbException("GbObject3D::addSurfaceTriangleSet - not implemented for "+(std::string)typeid(*this).name()); }

   virtual bool hasRaytracing() { return false; }
   virtual bool raytracingSupportsPointsInside() { return false; }
   //|r| must be 1! einheitsvector!!
   //return negativ value oder zero if no intersection
   virtual double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3) { throw UbException("GbObject3D::getIntersectionRaytraceFactor - not implemented"); }
#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<ObObject>(ar, *this);
   }
#endif //CAB_RCF
};
/*=========================================================================*/

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   SF_NO_CTOR(GbObject3D);
   UB_AUTO_RUN_NAMED(SF::registerType<GbObject3D>("GbObject3D") , SF_GbObject3D);
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived<ObObject, GbObject3D >() ), SF_GbObject3D_BD1 );
#endif //RCF_USE_SF_SERIALIZATION


#endif
