#ifndef DecreaseViscosityCoProcessor_H
#define DecreaseViscosityCoProcessor_H

#include <memory>

#include "CoProcessor.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"

#include "MuParser/include/muParser.h"


class DecreaseViscosityCoProcessor;
typedef std::shared_ptr<DecreaseViscosityCoProcessor> DecreaseViscosityCoProcessorPtr;

class UbScheduler;
class Grid3D;
class Communicator;

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
   DecreaseViscosityCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s,
      mu::Parser* nueFunc, std::shared_ptr<Communicator> comm);
   virtual ~DecreaseViscosityCoProcessor();
   //! calls collect PostprocessData.
   void process(double step); 
protected:
   //! resets the collision factor depending on the current timestep.
   void setViscosity(double step);  
   std::shared_ptr<Communicator>  comm;
private:
   mutable mu::value_type timeStep;
   mu::Parser* nueFunc;
};


#endif /* DecreaseViscosityCoProcessor_H_ */
