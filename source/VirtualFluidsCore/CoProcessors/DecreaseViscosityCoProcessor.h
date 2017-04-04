#ifndef DecreaseViscosityCoProcessor_H
#define DecreaseViscosityCoProcessor_H

#include "CoProcessor.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class DecreaseViscosityCoProcessor;
typedef boost::shared_ptr<DecreaseViscosityCoProcessor> DecreaseViscosityCoProcessorPtr;

//! \brief The class sets viscosity/collision factor according to a previously defined function in time. 
//! \details initialization in test case (example): 
//! \code{.cpp}
//! mu::Parser decrViscFunc;                       //define a mu-parser function 
//! decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");   //this function is time-dependent, the viscosity decreases a 1/t^2 
//! decrViscFunc.DefineConst("nue0", nueLB);       
//! decrViscFunc.DefineConst("c0", 0.1);           //constants such as c0 controll how fast the viscosity decreasis 
//! UbSchedulerPtr DecrViscSch(new UbScheduler()); //the CoProcessor is called according to a Scheduler
//! DecrViscSch->addSchedule(10,10,1000);          //in this case the viscosity is reset every 10 timesteps for the first 1000 timesteps 
//! DecreaseViscosityCoProcessor decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm); 
//! \endcode
//! \author Sonja Uphoff

class DecreaseViscosityCoProcessor: public CoProcessor 
{ 
public:
   DecreaseViscosityCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
      mu::Parser* nueFunc, CommunicatorPtr comm);
   virtual ~DecreaseViscosityCoProcessor();
   //! calls collect PostprocessData.
   void process(double step); 
protected:
   //! resets the collision factor depending on the current timestep.
   void setViscosity(double step);  
   CommunicatorPtr comm;
private:
   mutable mu::value_type timeStep;
   mu::Parser* nueFunc;
};


#endif /* DecreaseViscosityCoProcessor_H_ */
