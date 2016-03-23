//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef FEHALFDISC3D_H
#define FEHALFDISC3D_H

#include <vector>
#include <cmath>

#include <numerics/geometry3d/GbObject3D.h>
#include <basics/utilities/UbObserver.h>

class GbPoint3D;
class GbLine3D;
class GbTriangle3D;

class GbObject3DCreator;

class FeHalfDisc3D : public GbObject3D , public UbObserver 
{
public:
   FeHalfDisc3D();
	FeHalfDisc3D(const double& x1a,const double& x2a, const double& x3a, const double& x1b,const double& x2b, const double& x3b, const double& radius);
	FeHalfDisc3D(GbPoint3D* p1, GbPoint3D* p2, const double& radius);
	FeHalfDisc3D(GbLine3D* line, const double& rad);
	FeHalfDisc3D(FeHalfDisc3D* cylinder);
	~FeHalfDisc3D();    

	FeHalfDisc3D* clone() { return new FeHalfDisc3D(this); }
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

	double getX1Centroid();
	double getX1Minimum() ;
	double getX1Maximum() ;
	double getX2Centroid();
	double getX2Minimum() ;
	double getX2Maximum() ;
	double getX3Centroid();
	double getX3Minimum() ;
	double getX3Maximum() ;

	bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p); 
	bool isPointInGbObject3D(const double& x1p, const double& x2p, const double& x3p, bool& pointIsOnBoundary); 
   bool isCellInsideGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);
   bool isCellInsideOrCuttingGbObject3D(const double& x1a,const double& x2a,const double& x3a,const double& x1b,const double& x2b,const double& x3b);

	GbLine3D* createClippedLine3D(GbPoint3D& point1, GbPoint3D& point2);
   
   bool hasRaytracing() { return true; }
   /*|r| must be 1! einheitsvector!!*/
   double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3);

	std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   
	std::string toString();
	ObObjectCreator* getCreator();
	void write(UbFileOutput* out);
	void read(UbFileInput* in);

	//virtuelle Methoden von UbObserver
	void objectChanged(UbObservable* changedObject);
	void objectWillBeDeleted(UbObservable* objectForDeletion);

   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere

protected:
	void initCylinderType();

   GbLine3D* mLine;
	double    mRad;

	int cylinderType;

	//void berechneQuerschnittsWerte();
   static const int NOTPARALLELTOAXIS  = (1<<0); //1
   static const int X1PARALLEL         = (1<<1); //2
   static const int X2PARALLEL         = (1<<2); //4
   static const int X3PARALLEL         = (1<<3); //8
};

#endif   
