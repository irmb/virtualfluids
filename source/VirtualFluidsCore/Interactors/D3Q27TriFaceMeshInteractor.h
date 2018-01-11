#ifndef D3Q19AMRTRIFACEMESHINTERACTOR_H
#define D3Q19AMRTRIFACEMESHINTERACTOR_H

#include <string>
#include <vector>
#include <map>
#include <memory>

#include "D3Q27Interactor.h"
#include "CbArray3D.h"

class D3Q27TriFaceMeshInteractor;
typedef std::shared_ptr<D3Q27TriFaceMeshInteractor> D3Q27TriFaceMeshInteractorPtr;


class GbObject3D;
class Grid3D;
class BCAdapter;
class GbTriFaceMesh3D;
class Block3D;


class D3Q27TriFaceMeshInteractor : public D3Q27Interactor 
{
public:
   static const int STRESSNORMAL=0;
   static const int STRESSALTERNATIV=1;

   D3Q27TriFaceMeshInteractor();
   D3Q27TriFaceMeshInteractor(std::shared_ptr<Grid3D> grid, std::string name="D3Q27TriFaceMeshInteractor");
   D3Q27TriFaceMeshInteractor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, int type);
   D3Q27TriFaceMeshInteractor(std::shared_ptr<GbTriFaceMesh3D> triFaceMesh, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type);
   D3Q27TriFaceMeshInteractor(std::shared_ptr<GbTriFaceMesh3D> triFaceMesh, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter, int type, Interactor3D::Accuracy a);
   //D3Q27TriFaceMeshInteractor(GbTriFaceMesh3DPtr triFaceMesh, D3Q27BoundaryConditionAdapterPtr bcAdapter, int type, std::string name="D3Q27TriFaceMeshInteractor");

   ~D3Q27TriFaceMeshInteractor();

   virtual void initInteractor(const double& timeStep=0);
   virtual void initInteractor2(const double& timeStep=0);

   void updateInteractor(const double& timestep=0);

   void updateMovedGeometry(const double& timeStep=0);
   void setQs(const double& timeStep);
   void refineBlockGridToLevel(int level, double startDistance, double stopDistance);

   bool setDifferencesToGbObject3D(const std::shared_ptr<Block3D> block/*,const double& orgX1,const double& orgX2,const double& orgX3,const double& blockLengthX1,const double& blockLengthX2,const double& blockLengthX3, const double& timestep=0*/);

   void setRegardPointInObjectTest( bool opt ) { this->regardPIOTest = opt; }

   ObObject*        clone() { throw UbException(UB_EXARGS,"not implemented");	}
   ObObjectCreator* getCreator();

   UbTupleDouble3 getForces();
   UbTupleDouble3 getForcesTriangle();

   void setStressMode(int stressMode)                      { this->stressMode = stressMode;                         }
   void setUseHalfSpaceCheck(bool useHalfSpace )           { this->useHalfSpace = useHalfSpace;                     }
   //void setReinitWithStoredQs(bool reinitWithStoredQsFlag) { this->reinitWithStoredQsFlag = reinitWithStoredQsFlag; }

   void calculateForces();
   void calculateStresses(); 
   void calculateStressesAlternativ();            

   void calcStressesLine(UbTupleDouble6& stresses, const double& weight, const UbTupleDouble6& stvW, const UbTupleDouble6& stvE );
   void calcStressesFace(UbTupleDouble6& stresses, const double& weightX, const double& weightY, const UbTupleDouble6& stvSW, const UbTupleDouble6& stvSE, const UbTupleDouble6& stvNE, const UbTupleDouble6& stvNW );
   void calcStressesCube(UbTupleDouble6& stresses, const double& weightX, const double& weightY, const double& weightZ, const UbTupleDouble6& stvBSW, const UbTupleDouble6& stvBSE, const UbTupleDouble6& stvBNE, const UbTupleDouble6& stvBNW, const UbTupleDouble6& stvTSW, const UbTupleDouble6& stvTSE, const UbTupleDouble6& stvTNE, const UbTupleDouble6& stvTNW  );

   void calculatePressure(); 
   void calcPressureLine(double &p, const double& weight, const double& pW, const double& pE );
   void calcPressureFace(double &p, const double& weightX, const double& weightY, const double& pSW, const double& pSE, const double& pNE, const double& pNW );
   void calcPressureCube(double &p, const double& weightX, const double& weightY, const double& weightZ, const double& pBSW, const double& pBSE, const double& pBNE, const double& pBNW, const double& pTSW, const double& pTSE, const double& pTNE, const double& pTNW  );

   void   setForceShift(double forceshift)   { this->forceshift = forceshift; this->forceshiftpolicy = true; }
   void   setVelocityShift(double velocityshift)   { this->velocityshift = velocityshift; this->velocityshiftpolicy = true; }
   double getForceShift()     { return this->forceshift; }
   double getVelocityShift()  { return this->velocityshift; }
   bool   getForceShiftPolicy() { return forceshiftpolicy;}
   bool   getVelocityShiftPolicy() { return velocityshiftpolicy;}

   void clearBcNodeIndicesAndQsMap() { this->bcNodeIndicesAndQsMap.clear();}

   virtual std::string toString();


protected:
   int    stressMode;

   double forceshift;       
   double velocityshift;
   bool   forceshiftpolicy;
   bool   velocityshiftpolicy;
   bool   useHalfSpace;
   bool   regardPIOTest;

   void reinitWithStoredQs(const double& timeStep);
   //   bool reinitWithStoredQsFlag;
   std::map< std::shared_ptr<Block3D>, std::map < UbTupleInt3, std::vector< float > > > bcNodeIndicesAndQsMap;    //!!! es kann sein, dass in diesem interactor
   //an eine rpos eine BC gesetzt wurde, aber derselbe node in
   //in einem anderen in einen anderen Typ (z.B. Solid) geaendert
   //wurde --> es ist keine BC mehr an der stelle!

   enum SolidCheckMethod { ScanLine, PointInObject };

   enum FLAGS { BC_FLAG, UNDEF_FLAG, FLUID_FLAG, SOLID_FLAG, OLDSOLID_FLAG };
   void recursiveGridFill(CbArray3D<FLAGS>& flagfield, const short& xs, const short& ys, const short& zs, const FLAGS& type);
   void iterativeGridFill(CbArray3D<FLAGS>& flagfield, const short& xs, const short& ys, const short& zs, const FLAGS& type); 
};


#endif 
