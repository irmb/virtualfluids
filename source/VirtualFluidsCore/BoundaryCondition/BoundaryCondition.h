
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include <vector>
#include <string>

#include "D3Q27BoundaryCondition.h"
#include "D3Q27System.h"
#include "EsoTwist3D.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include <boost/shared_ptr.hpp>

class BoundaryCondition;
typedef boost::shared_ptr<BoundaryCondition> BoundaryConditionPtr;

class BoundaryCondition
{
public:
   enum Type
   {
      Velocity, Density, NoSlip, Slip 
   };
public:
   BoundaryCondition();
   virtual ~BoundaryCondition() {}
   
   void apply();
   virtual void addDistributions(DistributionArray3DPtr distributions) = 0;
   void addNode(int x1, int x2, int x3);
   void addBcPointer(D3Q27BoundaryConditionPtr bcPtr);
   void setCompressible(bool c);
   void setCollFactor(LBMReal cf);
   BoundaryCondition::Type getType();
   bool isPreCollision();
   virtual BoundaryConditionPtr clone()=0;
protected:
   virtual void applyBC() = 0;
   
   std::vector <int> nodeVector;
   std::vector <D3Q27BoundaryConditionPtr> bcVector;

   bool compressible;
   Type type;
   bool preCollision;

   D3Q27BoundaryConditionPtr bcPtr;
   DistributionArray3DPtr distributions;

   LBMReal collFactor;
   int x1, x2, x3;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   typedef LBMReal(*CalcFeqForDirFct)(const int& /*direction*/, const LBMReal& /*(d)rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   typedef  void(*CalcFeqFct)(LBMReal* const& /*feq/*[27]*/, const LBMReal& /*rho*/, const LBMReal& /*vx1*/, const LBMReal& /*vx2*/, const LBMReal& /*vx3*/);
   
   CalcFeqForDirFct calcFeqsForDirFct ;
   CalcMacrosFct    calcMacrosFct;
   CalcFeqFct       calcFeqFct;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & nodeVector;
      ar & bcVector;
      ar & compressible;
      ar & type;
      ar & distributions;
      ar & collFactor;
   }
};


#endif 
