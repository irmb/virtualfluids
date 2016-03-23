#ifndef D3Q27DensityBCAdapter_H
#define D3Q27DensityBCAdapter_H
        
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "basics/utilities/UbMath.h"
#include "basics/utilities/UbTuple.h"

#include "D3Q27BoundaryConditionAdapter.h"
#include "D3Q27BCFunction.h"

//*  D3Q27DensityBCAdapter                                                            */
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


class D3Q27DensityBCAdapter : public D3Q27BoundaryConditionAdapter
{
public:
   //constructors
   D3Q27DensityBCAdapter() { this->init(); }
   D3Q27DensityBCAdapter(const double& dens, const double& startTime=0.0, const double& endTime = D3Q27BCFunction::INFCONST );
   D3Q27DensityBCAdapter(const D3Q27BCFunction& densBC );
   D3Q27DensityBCAdapter(const std::vector< D3Q27BCFunction >& densBCs);
   D3Q27DensityBCAdapter(const mu::Parser& function, const double& startTime=0.0, const double& endTime = D3Q27BCFunction::INFCONST  );

   //------------- implements D3Q27BoundaryConditionAdapter ----- start
   std::string toString();
   ObObjectCreator* getCreator();

   void init(const D3Q27Interactor* const& interactor, const double& time=0);
   void update(const D3Q27Interactor* const& interactor, const double& time=0);

   void adaptBCForDirection( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 );
   void adaptBC( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 );

   double getDensity(const double& x1, const double& x2, const double& x3, const double& timeStep);

   //------------- implements D3Q27BoundaryConditionAdapter ----- end


protected:
   void init();
   
   //time dependency wird automatisch ueber D3Q27BCFunction Intervalle ermittelt!
   void setTimeDependent()   { (this->type |=   TIMEDEPENDENT);}
   void unsetTimeDependent() { (this->type &=  ~TIMEDEPENDENT);}
   
   void clear() { densBCs.clear(); }
   void setNodeDensity(const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& timestep);

private:
   mu::value_type x1, x2, x3; //brauch man nicht serialisieren!
   mu::value_type timeStep;   //brauch man nicht serialisieren!

   mu::Parser* tmpDensityFunction; //brauch man nicht serialisieren!
   
   std::vector<D3Q27BCFunction> densBCs;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<D3Q27BoundaryConditionAdapter>(*this);
      ar & densBCs;
   }
};

#endif 
