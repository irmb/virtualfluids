
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <memory>

#include "D3Q27System.h"

class DistributionArray3D;
class BCArray3D;
class BoundaryConditions;

class BCAlgorithm;
typedef std::shared_ptr<BCAlgorithm> BCAlgorithmPtr;

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
   static const char VelocityWithDensityBCAlgorithm = 7;
   static const char NonReflectingOutflowBCAlgorithm = 8;

public:
   BCAlgorithm();
   virtual ~BCAlgorithm() {}
   
   virtual void addDistributions(std::shared_ptr<DistributionArray3D> distributions) = 0;
   void setNodeIndex(int x1, int x2, int x3);
   void setBcPointer(std::shared_ptr<BoundaryConditions> bcPtr);
   void setCompressible(bool c);
   void setCollFactor(LBMReal cf);
   char getType();
   bool isPreCollision();
   virtual BCAlgorithmPtr clone() = 0;
   std::shared_ptr<BCArray3D> getBcArray();
   void setBcArray(std::shared_ptr<BCArray3D> bcarray);
   virtual void applyBC() = 0;

protected:
   bool compressible;
   char type;
   bool preCollision;

   std::shared_ptr<BoundaryConditions> bcPtr;
   std::shared_ptr<DistributionArray3D> distributions;
   std::shared_ptr<BCArray3D> bcArray;

   LBMReal collFactor;
   int x1, x2, x3;

   LBMReal compressibleFactor;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   typedef LBMReal(*CalcFeqForDirFct)(const int& /*direction*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   typedef  void(*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/, const LBMReal& /*rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   
   CalcFeqForDirFct calcFeqsForDirFct ;
   CalcMacrosFct    calcMacrosFct;
   CalcFeqFct       calcFeqFct; 
};


#endif 
