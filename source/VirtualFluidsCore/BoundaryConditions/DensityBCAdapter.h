#ifndef DensityBCAdapter_H
#define DensityBCAdapter_H
        
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "basics/utilities/UbMath.h"
#include "basics/utilities/UbTuple.h"

#include "BCAdapter.h"
#include "BCFunction.h"

//*  DensityBCAdapter                                                            */
//*                                                                         */
//**
//<BR><BR>
//@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
//@version 1.0 - 06.09.06
//*/ 
//
//*
//usage: ...
//*/


class DensityBCAdapter : public BCAdapter
{
public:
   //constructors
   DensityBCAdapter() { this->init(); }
   DensityBCAdapter(const double& dens, const double& startTime=0.0, const double& endTime = BCFunction::INFCONST );
   DensityBCAdapter(const BCFunction& densBC );
   DensityBCAdapter(const std::vector< BCFunction >& densBCs);
   DensityBCAdapter(const mu::Parser& function, const double& startTime=0.0, const double& endTime = BCFunction::INFCONST  );

   //------------- implements D3Q27BoundaryConditionAdapter ----- start
   std::string toString();
   ObObjectCreator* getCreator();

   void init(const D3Q27Interactor* const& interactor, const double& time=0);
   void update(const D3Q27Interactor* const& interactor, const double& time=0);

   void adaptBCForDirection( const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 );
   void adaptBC( const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 );

   double getDensity(const double& x1, const double& x2, const double& x3, const double& timeStep);

   //------------- implements D3Q27BoundaryConditionAdapter ----- end


protected:
   void init();
   
   //time dependency wird automatisch ueber D3Q27BCFunction Intervalle ermittelt!
   void setTimeDependent()   { (this->type |=   TIMEDEPENDENT);}
   void unsetTimeDependent() { (this->type &=  ~TIMEDEPENDENT);}
   
   void clear() { densBCs.clear(); }
   void setNodeDensity(const D3Q27Interactor& interactor, BoundaryConditionsPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& timestep);

private:
   mu::value_type x1, x2, x3; //brauch man nicht serialisieren!
   mu::value_type timeStep;   //brauch man nicht serialisieren!

   mu::Parser* tmpDensityFunction; //brauch man nicht serialisieren!
   
   std::vector<BCFunction> densBCs;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<BCAdapter>(*this);
      ar & densBCs;
   }
};

#endif 
