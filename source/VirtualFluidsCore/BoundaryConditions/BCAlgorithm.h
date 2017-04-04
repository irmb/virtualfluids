
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <vector>
#include <string>

#include "BoundaryConditions.h"
#include "D3Q27System.h"
#include "EsoTwist3D.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include <boost/shared_ptr.hpp>

class BCAlgorithm;
typedef boost::shared_ptr<BCAlgorithm> BCAlgorithmPtr;

class BCAlgorithm
{
public:
   static const char VelocityBCAlgorithm = 0;
   static const char EqDensityBCAlgorithm = 1;
   static const char NonEqDensityBCAlgorithm = 2;
   static const char NoSlipBCAlgorithm = 3;
   static const char SlipBCAlgorithm = 4;
   static const char HighViscosityNoSlipBCAlgorithm = 5;
   static const char ThinWallNoSlipBCAlgorithm = 6;
   static const char NonReflectingVelocityBCAlgorithm = 7;
   static const char NonReflectingDensityBCAlgorithm = 8;
public:
   BCAlgorithm();
   virtual ~BCAlgorithm() {}
   
   void apply();
   virtual void addDistributions(DistributionArray3DPtr distributions) = 0;
   void addNode(int x1, int x2, int x3);
   void addBcPointer(BoundaryConditionsPtr bcPtr);
   void setCompressible(bool c);
   void setCollFactor(LBMReal cf);
   char getType();
   bool isPreCollision();
   virtual BCAlgorithmPtr clone()=0;
   void clearData();
protected:
   virtual void applyBC() = 0;
   
   std::vector <int> nodeVector;
   std::vector <BoundaryConditionsPtr> bcVector;

   bool compressible;
   char type;
   bool preCollision;

   BoundaryConditionsPtr bcPtr;
   DistributionArray3DPtr distributions;

   LBMReal collFactor;
   int x1, x2, x3;

   LBMReal compressibleFactor;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   typedef LBMReal(*CalcFeqForDirFct)(const int& /*direction*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   typedef  void(*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/, const LBMReal& /*rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   
   CalcFeqForDirFct calcFeqsForDirFct ;
   CalcMacrosFct    calcMacrosFct;
   CalcFeqFct       calcFeqFct;

private:
   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & nodeVector;
   //   ar & bcVector;
   //   ar & compressible;
   //   ar & type;
   //   ar & distributions;
   //   ar & collFactor;
   //}
};


#endif 
