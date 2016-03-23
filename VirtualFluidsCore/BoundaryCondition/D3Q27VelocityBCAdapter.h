//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef D3Q27VELOCITYADAPTER_H
#define D3Q27VELOCITYADAPTER_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif

#include <basics/utilities/UbInfinity.h>
#include <basics/utilities/UbFileOutput.h>
#include <basics/utilities/UbFileInput.h>

class UbFileOutput;
class UbFileInput;

#include <D3Q27BoundaryConditionAdapter.h>
#include <D3Q27BCFunction.h>

//example:
//        vector<D3Q27BCFunction> vx1BCs,vx2BCs,vx3BCs;
//        vx1BCs.push_back(D3Q27BCFunction(0.01 , 0  , 100) );   //t=[0  ..100[ -> vx1 = 0.01
//        vx1BCs.push_back(D3Q27BCFunction(0.004, 100, 200) );   //t=[100..200[ -> vx1 = 0.004
//        vx1BCs.push_back(D3Q27BCFunction(0.03 , 200, 400) );   //t=[200..400] -> vx1 = 0.03
// 
//        vx2BCs.push_back(D3Q27BCFunction(0.02 , 0  , 200) );   //t=[0  ..200[ -> vx2 = 0.02
//        vx2BCs.push_back(D3Q27BCFunction(0.002, 200, 300) );   //t=[200..300[ -> vx2 = 0.002
//        vx2BCs.push_back(D3Q27BCFunction(0.043, 300, 600) );   //t=[300..600] -> vx2 = 0.043
//        
//        D3Q27VelocityBCAdapter bcAdapter(vx1BCs,vx2BCs,vx3BCs);
//        bcAdapter.setTimePeriodic(); //->  t=[0  ..100[ -> vx1 = 0.01
//                                     //    t=[100..200[ -> vx1 = 0.004
//                                     //    t=[200..400[ -> vx1 = 0.03
//                                     //    t=[400..500[ -> vx1 = 0.01
//                                     //    t=[500..600[ -> vx1 = 0.004
//                                     //    t=[600..800[ -> vx1 = 0.03  ...
//                                     //    t=[0  ..200[ -> vx2 = 0.02
//                                     //    t=[200..300[ -> vx2 = 0.002
//                                     //    t=[300..600] -> vx2 = 0.043
//                                     //    t=[600..800[ -> vx2 = 0.02
//                                     //    t=[800..900[ -> vx2 = 0.002
//                                     //    t=[900..1200]-> vx2 = 0.043  ...
//
// example parabolic inflow:
//    mu::Parser fct;
//    fct.SetExpr("max(vmax*(1.0-4.0*((x2-x2_vmax)^2+(x3-x3_vmax)^2)/H^2),0.0)"); //paraboloid (mit vmax bei (0/x2_vmax/x3_vmax) 
//    fct.DefineConst("x2Vmax", 0.0            ); //x2-Pos für vmax
//    fct.DefineConst("x3Vmax", 0.0            ); //x3-Pos für vmax
//    fct.DefineConst("H"     , rohrDurchmesser);
//    fct.DefineConst("vmax"  , vmax           );
//    D3Q27VelocityBCAdapter velBC(true, false ,false ,fct, 0, D3Q27BCFunction::INFCONST);

/*=========================================================================*/
/*  D3Q27VelocityBCAdapter                                                 */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 06.09.06
*/ 

class D3Q27VelocityBCAdapter : public D3Q27BoundaryConditionAdapter
{
public:
   //constructors
   D3Q27VelocityBCAdapter() { this->init(); }
   
   D3Q27VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const D3Q27BCFunction& velVxBC );

   D3Q27VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const mu::Parser& function, const double& startTime, const double& endTime  );

   D3Q27VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const mu::Parser& function1, const mu::Parser& function2, const mu::Parser& function3, const double& startTime, const double& endTime );
   
   D3Q27VelocityBCAdapter(const bool& vx1, const bool& vx2, const bool& vx3, const std::string& functionstring, const double& startTime, const double& endTime );

   D3Q27VelocityBCAdapter(const D3Q27BCFunction& velBC, bool x1Dir, bool x2Dir, bool x3Dir);

   D3Q27VelocityBCAdapter(const D3Q27BCFunction& velVx1BC, const D3Q27BCFunction& velVx2BC, const D3Q27BCFunction& velVx3BC);

   D3Q27VelocityBCAdapter(const std::vector< D3Q27BCFunction >& velVx1BCs, const std::vector< D3Q27BCFunction >& velVx2BCs, const std::vector< D3Q27BCFunction >& velVx3BCs);

   D3Q27VelocityBCAdapter(const double& vx1, const double& vx1StartTime, const double& vx1EndTime,
                          const double& vx2, const double& vx2StartTime, const double& vx2EndTime,
                          const double& vx3, const double& vx3StartTime, const double& vx3EndTime);

   D3Q27VelocityBCAdapter(const std::string& vx1Function, const double& vx1StartTime, const double& vx1EndTime,
                          const std::string& vx2Function, const double& vx2StartTime, const double& vx2EndTime,
                          const std::string& vx3Function, const double& vx3StartTime, const double& vx3EndTime ); 

   //methods
   void setTimePeriodic()    { (this->type |=   TIMEPERIODIC); }
   void unsetTimePeriodic()  { (this->type &=  ~TIMEPERIODIC); }
   bool isTimePeriodic()     { return ((this->type & TIMEPERIODIC) ==  TIMEPERIODIC); }

   //folgendes ist fuer moving objects gedadacht... 
   void setNewVelocities(const double& vx1, const double& vx1StartTime, const double& vx1EndTime,
                         const double& vx2, const double& vx2StartTime, const double& vx2EndTime,
                         const double& vx3, const double& vx3StartTime, const double& vx3EndTime);

      
   //------------- implements D3Q27BoundaryConditionAdapter ----- start
   std::string toString();
   
   void init(const D3Q27Interactor* const& interactor, const double& time=0);
   void update(const D3Q27Interactor* const& interactor, const double& time=0);

   void adaptBCForDirection( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time=0 );
   void adaptBC( const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time=0 );

   //------------- implements D3Q27BoundaryConditionAdapter ----- end

   UbTupleDouble3 getVelocity(const double& x1, const double& x2, const double& x3, const double& timeStep) const;


protected:
   void init();
   void init(std::vector<D3Q27BCFunction>& vxBCs);

   //time dependency wird automatisch ueber D3Q27BCFunction Intervalle ermittelt!
   void setTimeDependent()   { (this->type |=   TIMEDEPENDENT); }
   void unsetTimeDependent() { (this->type &=  ~TIMEDEPENDENT); }

   void clear() { vx1BCs.clear(); vx2BCs.clear();  vx3BCs.clear(); this->init(); }
   void setNodeVelocity(const D3Q27Interactor& interactor, D3Q27BoundaryConditionPtr bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& timestep);

private:
   mutable mu::value_type x1, x2, x3;
   mutable mu::value_type timeStep;

   mu::Parser* tmpVx1Function;
   mu::Parser* tmpVx2Function;
   mu::Parser* tmpVx3Function;

   std::vector<D3Q27BCFunction> vx1BCs;
   std::vector<D3Q27BCFunction> vx2BCs;
   std::vector<D3Q27BCFunction> vx3BCs;

private:
   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<D3Q27BoundaryConditionAdapter>(*this);
   }
};

#endif
