//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBCYLINDER3D_H
#define GBCYLINDER3D_H

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <vector>
#include <cmath>

#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/GbLine3D.h>
#include <basics/utilities/UbObserver.h>

class GbPoint3D;
class GbLine3D;
class GbTriangle3D;

class GbObject3DCreator;

#include <PointerDefinitions.h>
class GbCylinder3D;
typedef SPtr<GbCylinder3D> GbCylinder3DPtr;


class GbCylinder3D : public GbObject3D , public UbObserver 
{
public:
   GbCylinder3D();
	GbCylinder3D(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b, const double& radius);
	GbCylinder3D(GbPoint3D* p1, GbPoint3D* p2, const double& radius);
	GbCylinder3D(GbLine3D* line, const double& rad);
	GbCylinder3D(GbCylinder3D* cylinder);
	~GbCylinder3D();    

	GbCylinder3D* clone() { return new GbCylinder3D(this); }
	void finalize();

	double     getRadius() { return this->mRad; };
	GbLine3D*  getLine() {return mLine;}
	GbPoint3D* getPoint1();
	GbPoint3D* getPoint2();

	void setRadius(const double& radius);
	void setLine(GbLine3D* line);
	void setPoint1(const double& x1, const double& x2, const double& x3);
	void setPoint2(const double& x1, const double& x2, const double& x3);

	bool isParallelToX1Axis() { return((this->cylinderType & X1PARALLEL        )    ==  X1PARALLEL        );}
	bool isParallelToX2Axis() { return((this->cylinderType & X2PARALLEL        )    ==  X2PARALLEL        );}
	bool isParallelToX3Axis() { return((this->cylinderType & X3PARALLEL        )    ==  X3PARALLEL        );}
	bool isNotParallelToAxis(){ return((this->cylinderType & NOTPARALLELTOAXIS )    ==  NOTPARALLELTOAXIS );}

	double getHeight(); 

	void scale(const double& sx1, const double& sx2, const double& sx3);

   void translate(const double& x1, const double& x2, const double& x3) 
   {
      this->mLine->translate( x1, x2, x3 );
      this->calculateValues();
      //this->notifyObserversObjectChanged();
   }

   double getX1Centroid() { return centerX1; }
   double getX1Minimum()  { return minX1;    }
	double getX1Maximum()  { return maxX1;    }
	double getX2Centroid() { return centerX2; }
	double getX2Minimum()  { return minX2;    }
	double getX2Maximum()  { return maxX2;    }
	double getX3Centroid() { return centerX3; }
	double getX3Minimum()  { return minX3;    }
	double getX3Maximum()  { return maxX3;    }

	bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p); 
	bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary); 
   bool isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);

	GbLine3D* createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2);
   
   //SG ausdokumentieren, da der nur unendlcihe Zylinder macht ...
   //bool hasRaytracing() { return true; }
   bool hasRaytracing() { return false; }
   bool raytracingSupportsPointsInside() { return true; }
   
   
   /*|r| must be 1! einheitsvector!!*/
   double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3);

	std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   void addSurfaceTriangleSet(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles);
   void addSurfaceTriangleSetSegments(std::vector<UbTupleFloat3>& nodes, std::vector<UbTupleInt3>& triangles, int segmentsRound, int segmentsHeight );

	std::string toString();
	ObObjectCreator* getCreator();
	void write(UbFileOutput* out);
	void read(UbFileInput* in);

	//virtuelle Methoden von UbObserver
	void objectChanged(UbObservable* changedObject);
	void objectWillBeDeleted(UbObservable* objectForDeletion);

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & mLine;
      ar & mRad;
      ar & cylinderType;
      
      if( ArchiveTools::isReading(ar) )
         this->calculateValues();
   }
#endif //CAB_RCF
   
   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere


protected:
   void calculateValues();

   GbLine3D* mLine;
	double    mRad;

   double minX1, minX2, minX3;
   double maxX1, maxX2, maxX3;
   double centerX1, centerX2, centerX3;

	int cylinderType;

	//void berechneQuerschnittsWerte();
   static const int NOTPARALLELTOAXIS  = (1<<0); //1
   static const int X1PARALLEL         = (1<<1); //2
   static const int X2PARALLEL         = (1<<2); //4
   static const int X3PARALLEL         = (1<<3); //8
};

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   UB_AUTO_RUN_NAMED(   SF::registerType<GbCylinder3D >("GbCylinder3D")           , SF_GbCylinder3D     );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbCylinder3D >() ), SF_GbCylinder3D_BD1 );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< UbObserver, GbCylinder3D>()  ), SF_GbCylinder3D_BD2 );
#endif //RCF_USE_SF_SERIALIZATION

#endif   
