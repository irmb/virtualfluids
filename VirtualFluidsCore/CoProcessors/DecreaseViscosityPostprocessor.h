#ifndef DecreaseViscosityPOSTPROCESSOR_H
#define DecreaseViscosityPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class DecreaseViscosityPostprocessor;
typedef boost::shared_ptr<DecreaseViscosityPostprocessor> DecreaseViscosityPostprocessorPtr;

//! \brief The class sets viscosity/collision factor according to a previously defined function in time. 
//! \details initialization in test case (example): 
//! \code{.cpp}
//! mu::Parser decrViscFunc;                       //define a mu-parser function 
//! decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");   //this function is time-dependent, the viscosity decreases a 1/t^2 
//! decrViscFunc.DefineConst("nue0", nueLB);       
//! decrViscFunc.DefineConst("c0", 0.1);           //constants such as c0 controll how fast the viscosity decreasis 
//! UbSchedulerPtr DecrViscSch(new UbScheduler()); //the postprocessor is called according to a Scheduler
//! DecrViscSch->addSchedule(10,10,1000);          //in this case the viscosity is reset every 10 timesteps for the first 1000 timesteps 
//! DecreaseViscosityPostprocessor decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm); 
//! \endcode
//! \author Sonja Uphoff

class DecreaseViscosityPostprocessor: public Postprocessor 
{ 
public:
   DecreaseViscosityPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
      mu::Parser* nueFunc, CommunicatorPtr comm);
   virtual ~DecreaseViscosityPostprocessor();
   //! calls collect PostprocessData.
   void update(double step); 
protected:
   //! resets the collision factor depending on the current timestep.
   void collectPostprocessData(double step);  
   CommunicatorPtr comm;
private:
   mutable mu::value_type timeStep;
   mu::Parser* nueFunc;
};


#endif /* DecreaseViscosityPOSTPROCESSOR_H_ */
