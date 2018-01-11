//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef D3Q27INTERACTOR_H
#define D3Q27INTERACTOR_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>

#include "UbException.h"
#include "UbTuple.h"
#include "GbPoint3D.h"
#include "Interactor3D.h"
#include "D3Q27System.h"

class D3Q27Interactor;
typedef std::shared_ptr<D3Q27Interactor> D3Q27InteractorPtr;

class BCAdapter;
class Block3D;
class Grid3D;
class GbObject3D;

typedef std::map<std::shared_ptr<Block3D>, std::set< std::vector<int> > > BcNodeIndicesMap;
typedef std::map<std::shared_ptr<Block3D>, std::set< UbTupleInt3 > > SolidNodeIndicesMap;

class D3Q27Interactor : public Interactor3D 
{
public:
   D3Q27Interactor();
   D3Q27Interactor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, int type);
   D3Q27Interactor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter,  int type);
   D3Q27Interactor(std::shared_ptr<GbObject3D> geoObject3D, std::shared_ptr<Grid3D> grid, std::shared_ptr<BCAdapter> bcAdapter,  int type, Interactor3D::Accuracy a);

   virtual ~D3Q27Interactor();

   void setRelevantForForces(const bool& value) {  this->relevantForForces = value; }
   bool isRelevantForForces() { return this->relevantForForces; }

   virtual void addBCAdapter(const std::shared_ptr<BCAdapter> bcAdapter) { bcAdapters.push_back(bcAdapter); }
   void deleteBCAdapter() { bcAdapters.clear(); }

 
   virtual void initInteractor(const double& timeStep=0);
   void updateInteractor(const double& timestep=0); 

   void setReinitWithStoredQs(bool reinitWithStoredQsFlag) { this->reinitWithStoredQsFlag = reinitWithStoredQsFlag; }
   
   void removeSolidBlocks() { Interactor3D::removeSolidBlocks(); solidNodeIndicesMap.clear(); }
   void removeBcBlocks() { Interactor3D::removeBcBlocks(); bcNodeIndicesMap.clear(); }

   bool setDifferencesToGbObject3D(const std::shared_ptr<Block3D> block/*, const double& x1, const double& x2, const double& x3, const double& blockLengthX1, const double& blockLengthX2, const double& blockLengthX3, const double& timestep=0*/);

   ObObject* clone() { throw UbException(UB_EXARGS,"not implemented");	}
   ObObjectCreator* getCreator();


   void writeValidationAVSFile(std::string filename);  
   virtual std::vector< std::pair<GbPoint3D,GbPoint3D> >  getQsLineSet();

   void addQsLineSet(std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines);

   const BcNodeIndicesMap& getBcNodeIndicesMap() const { return bcNodeIndicesMap; }

protected:
   bool relevantForForces;
   bool reinitWithStoredQsFlag;

   std::vector<std::shared_ptr<BCAdapter> > bcAdapters;
   

   SolidNodeIndicesMap solidNodeIndicesMap;
   BcNodeIndicesMap bcNodeIndicesMap;
                                                                         //!!! es kann sein, dass in diesem interactor
                                                                         //an eine rpos eine BC gesetzt wurde, aber derselbe node in
                                                                         //in einem anderen in einen anderen Typ (z.B. Solid) geaendert
                                                                         //wurde --> es ist keine BC mehr an der stelle!
   
   void   initRayVectors();
   double rayX1[D3Q27System::FENDDIR+1];
   double rayX2[D3Q27System::FENDDIR+1];
   double rayX3[D3Q27System::FENDDIR+1];

};


#endif
