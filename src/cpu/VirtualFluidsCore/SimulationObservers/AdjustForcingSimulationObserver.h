#ifndef D3Q27ADJUSTFORCINGSimulationObserver_H
#define D3Q27ADJUSTFORCINGSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"
#include "lbm/constants/D3Q27.h"

namespace vf::parallel {class Communicator;}
class UbScheduler;
class Grid3D;
class IntegrateValuesHelper;

//! \brief   Computes forcing such that a given velocity (vx1Targed) is reached inside an averaging domain (h1).
//! \details Algorithm based on PID controller (proportional_integral_derivative controller). The parameters of PID
//! controller estimation based on Ziegler�Nichols method.
//!          Integrate values helper, scheduler must be set in test case.
//! \author: Konstantin Kutscher

class AdjustForcingSimulationObserver : public SimulationObserver
{
public:
    AdjustForcingSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                             SPtr<IntegrateValuesHelper> integrateValues, real vTarged, std::shared_ptr<vf::parallel::Communicator> comm);
    //!< calls collect PostprocessData
    void update(real step) override;

protected:
    //!< object that can compute spacial average values in 3D-subdomain.
    SPtr<IntegrateValuesHelper> integrateValues;
    //!< compares velocity in integrateValues with target velocity and adjusts forcing accordingly.
    void collectData(real step);
    std::shared_ptr<vf::parallel::Communicator> comm;

private:
    real vx1Targed; //!< target velocity.
    real forcing;   //!< forcing at previous update step.
    real cellsVolume;
    real vx1Average;
    bool root;
    real Kpcrit; // Kp critical
    real Tcrit;  // the oscillation period
    real Tn;
    real Tv;
    real e;
    real Ta;
    real Kp;
    real Ki;
    real Kd;
    real y;
    real esum;
    real eold;
    // std::vector<CalcNodes> cnodes;
    std::string path;
};

#endif /* D3Q27RHODIFFERENCESimulationObserver_H_ */
