//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef D3Q27INTERACTOR_H
#define D3Q27INTERACTOR_H

#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <cmath>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/set.hpp>

#include <boost/shared_ptr.hpp>

class D3Q27Interactor;
typedef boost::shared_ptr<D3Q27Interactor> D3Q27InteractorPtr;


#include "UbException.h"
#include "UbTuple.h"
#include "ObFactory.h"
#include "CoordinateTransformation3D.h"
#include "GbPoint3D.h"
#include "Interactor3D.h"
#include "BCArray3D.h"
#include "BCAdapter.h"
#include "BoundaryConditions.h"
#include "D3Q27System.h"

class UbFileInput;
class UbFileOutput;
class GbObject3D;

//////////////////////////////////////////////////////////////////////////
class D3Q27Interactor : public Interactor3D 
{
public:
   D3Q27Interactor();
   D3Q27Interactor(GbObject3DPtr geoObject3D, Grid3DPtr grid, int type);
   D3Q27Interactor(GbObject3DPtr geoObject3D, Grid3DPtr grid, BCAdapterPtr bcAdapter,  int type);
   D3Q27Interactor(GbObject3DPtr geoObject3D, Grid3DPtr grid, BCAdapterPtr bcAdapter,  int type, Interactor3D::Accuracy a);

   virtual ~D3Q27Interactor();

   void setRelevantForForces(const bool& value) {  this->relevantForForces = value; }
   bool isRelevantForForces() { return this->relevantForForces; }
   //UbTupleDouble3 getForces();
   //UbTupleDouble3 getForces(Patch3DPtr patch);

   virtual void addBCAdapter(const BCAdapterPtr bcAdapter) { bcAdapterVector.push_back(bcAdapter); }
   void deleteBCAdapter() { bcAdapterVector.clear(); }
   //virtual std::vector< MbSmartPtr<D3Q27BoundaryConditionAdapter> > getBcAdapters() { return bcAdapterVector; }

 
   virtual void initInteractor(const double& timeStep=0);
   void updateInteractor(const double& timestep=0); 
   //virtual void updateMovedGeometry(const double& timeStep=0); 
   //void updateNewNodes();
   void setReinitWithStoredQs(bool reinitWithStoredQsFlag) { this->reinitWithStoredQsFlag = reinitWithStoredQsFlag; }
   
   void removeSolidBlocks() { Interactor3D::removeSolidBlocks(); solidNodeIndicesMap.clear(); }
   void removeTransBlocks() { Interactor3D::removeTransBlocks(); transNodeIndicesMap.clear(); }
   virtual void removeBoundaryInformationOnTransNodes();

   bool setDifferencesToGbObject3D(const Block3DPtr block/*, const double& x1, const double& x2, const double& x3, const double& blockLengthX1, const double& blockLengthX2, const double& blockLengthX3, const double& timestep=0*/);

   ObObject* clone() { throw UbException(UB_EXARGS,"not implemented");	}
   ObObjectCreator* getCreator();

   //std::string toString();

   //------------- implements CAB serialization ----- start
   //void write(UbFileOutput* out);
   //void read(UbFileInput* in);
   //std::string writeTransNodes(std::string filename);
   //std::string writeTransNodesAsTriangles(std::string filename);
   void writeValidationAVSFile(std::string filename);  
   virtual std::vector< std::pair<GbPoint3D,GbPoint3D> >  getQsLineSet();

   void addQsLineSet(std::vector<UbTupleFloat3 >& nodes, std::vector<UbTupleInt2 >& lines);

   const std::map<Block3DPtr, std::set< std::vector<int> > > & getTransNodeIndicesMap() { return transNodeIndicesMap; }

protected:
   bool relevantForForces;
   bool reinitWithStoredQsFlag;

   //std::vector< MbSmartPtr<D3Q27BoundaryConditionAdapter> > bcAdapterVector;
   std::vector<BCAdapterPtr> bcAdapterVector;
   
   typedef UbTuple<int,int,int,long long>  UbTupleInt3LongLong;

   std::vector<Block3D*> oldSolidBlockSet;
   std::map<Block3DPtr, std::set< UbTupleInt3 > > oldSolidNodeIndicesMap;
   std::map<Block3DPtr, std::set< UbTupleInt3LongLong > > oldTransNodeIndicesMap;

   std::map<Block3DPtr, std::set< UbTupleInt3 > > solidNodeIndicesMap;  

   std::map<Block3DPtr, std::set< std::vector<int> > > transNodeIndicesMap;
   //std::map<Block3DPtr, std::set< UbTupleInt3 > > transNodeIndicesMap;//!!! es kann sein, dass in diesem interactor
                                                                         //an eine rpos eine BC gesetzt wurde, aber derselbe node in
                                                                         //in einem anderen in einen anderen Typ (z.B. Solid) geaendert
                                                                         //wurde --> es ist keine BC mehr an der stelle!
   
   void   initRayVectors();
   double rayX1[D3Q27System::FENDDIR+1];
   double rayX2[D3Q27System::FENDDIR+1];
   double rayX3[D3Q27System::FENDDIR+1];

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<Interactor3D>(*this);
      ar & transNodeIndicesMap;
      ar & bcAdapterVector;
   }
};


#endif
