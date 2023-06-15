#ifndef DecreaseViscositySimulationObserver_H
#define DecreaseViscositySimulationObserver_H

#include <PointerDefinitions.h>

#include "SimulationObserver.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"

#include "muParser.h"

class UbScheduler;
class Grid3D;
namespace vf::mpi {class Communicator;}

//! \brief The class sets viscosity/collision factor according to a previously defined function in time.
//! \details initialization in test case (example):
//! \code{.cpp}
//! mu::Parser decrViscFunc;                       //define a mu-parser function
//! decrViscFunc.SetExpr("nue0+c0/(t+1)/(t+1)");   //this function is time-dependent, the viscosity decreases a 1/t^2
//! decrViscFunc.DefineConst("nue0", nueLB);
//! decrViscFunc.DefineConst("c0", 0.1);           //constants such as c0 controll how fast the viscosity decreasis
//! SPtr<UbScheduler> DecrViscSch(new UbScheduler()); //the SimulationObserver is called according to a Scheduler
//! DecrViscSch->addSchedule(10,10,1000);          //in this case the viscosity is reset every 10 timesteps for the
//! first 1000 timesteps DecreaseViscositySimulationObserver decrViscPPPtr(grid, DecrViscSch,&decrViscFunc, comm); \endcode
//! \author Sonja Uphoff

class DecreaseViscositySimulationObserver : public SimulationObserver
{
public:
    DecreaseViscositySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, mu::Parser *nueFunc, std::shared_ptr<vf::mpi::Communicator> comm);
    ~DecreaseViscositySimulationObserver() override;
    //! calls collect PostprocessData.
    void update(real step) override;

protected:
    //! resets the collision factor depending on the current timestep.
    void setViscosity(real step);
    std::shared_ptr<vf::mpi::Communicator> comm;

private:
    mutable mu::value_type timeStep;
    mu::Parser *nueFunc;
};

#endif /* DecreaseViscositySimulationObserver_H_ */
