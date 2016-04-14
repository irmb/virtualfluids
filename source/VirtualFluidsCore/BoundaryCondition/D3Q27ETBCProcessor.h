#ifndef D3Q27ETBCPROCESSSOR_H
#define D3Q27ETBCPROCESSSOR_H

#include "BCProcessor.h"
#include "D3Q27BoundaryCondition.h"
#include "EsoTwist3D.h"
#include "BCArray3D.h"
#include "basics/container/CbArray4D.h"
#include "basics/container/CbArray3D.h"
//#include "BoundaryCondition.h"
#include <vector>

#include <boost/serialization/base_object.hpp>

class D3Q27ETBCProcessor;
typedef boost::shared_ptr<D3Q27ETBCProcessor> D3Q27ETBCProcessorPtr;

class D3Q27ETBCProcessor : public BCProcessor
{
public:
   D3Q27ETBCProcessor();
   D3Q27ETBCProcessor(LBMKernel3DPtr kernel);
   virtual ~D3Q27ETBCProcessor();
   //virtual void applyBC();
   virtual BCArray3D<D3Q27BoundaryCondition>& getBCArray();
   virtual BCProcessorPtr clone(LBMKernel3DPtr kernel);

   void addBC(BoundaryConditionPtr bc);
   BoundaryConditionPtr getBC(BoundaryCondition::Type type);
   void applyPreCollisionBC();
   void applyPostCollisionBC();
   //void init();
protected:
   std::vector<BoundaryConditionPtr> preBC;
   std::vector<BoundaryConditionPtr> postBC;

   //int minX1, minX2, minX3, maxX1, maxX2, maxX3;
   //D3Q27BoundaryConditionPtr bcPtr;
   //LBMReal f[D3Q27System::ENDF+1];
   //LBMReal ftemp[D3Q27System::ENDF+1];
   //LBMReal feq[D3Q27System::ENDF+1];
   //LBMReal rho, vx1, vx2, vx3, rho0;

   //void init(LBMKernel3DPtr kernel);
   BCArray3D<D3Q27BoundaryCondition> bcArray;
   //EsoTwist3DPtr distributions;
   //int ghostLayerWitdh;
   //LBMReal collFactor;
   //bool compressible;
   //EsoTwist3DPtr distributionsTemp;

   //std::vector <int> sizeVector;
   //std::vector <int> nodeVector;
   //std::vector <D3Q27BoundaryConditionPtr> bcVector;

   //typedef void(*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   //typedef LBMReal(*CalcFeqForDirFct)(const int& /*direction*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   //typedef  void(*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/, const LBMReal& /*rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   //CalcFeqForDirFct calcFeqsForDirFct;
   //CalcMacrosFct    calcMacrosFct;
   //CalcFeqFct       calcFeqFct;
private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BCProcessor>(*this);
      ar & bcArray;
      //ar & distributions;
      //ar & ghostLayerWitdh;
      //ar & collFactor;
      //ar & compressible;
      //ar & minX1;
      //ar & minX2;
      //ar & minX3;
      //ar & maxX1;
      //ar & maxX2;
      //ar & maxX3;
      //ar & distributionsTemp;
      ar & preBC;
      ar & postBC;
   }
};

#endif
