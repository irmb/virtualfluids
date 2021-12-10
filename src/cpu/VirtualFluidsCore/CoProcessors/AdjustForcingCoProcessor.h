#ifndef D3Q27ADJUSTFORCINGCoProcessor_H
#define D3Q27ADJUSTFORCINGCoProcessor_H

#include <PointerDefinitions.h>
#include <string>

#include "CoProcessor.h"

namespace vf::mpi {class Communicator;}
class UbScheduler;
class Grid3D;
class IntegrateValuesHelper;

//! \brief   Computes forcing such that a given velocity (vx1Targed) is reached inside an averaging domain (h1).
//! \details Algorithm based on PID controller (proportional_integral_derivative controller). The parameters of PID
//! controller estimation based on Ziegler�Nichols method.
//!          Integrate values helper, scheduler must be set in test case.
//! \author: Konstantin Kutscher

class AdjustForcingCoProcessor : public CoProcessor
{
public:
    AdjustForcingCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                             SPtr<IntegrateValuesHelper> integrateValues, double vTarged, std::shared_ptr<vf::mpi::Communicator> comm);
    //!< calls collect PostprocessData
    void process(double step) override;

protected:
    //!< object that can compute spacial average values in 3D-subdomain.
    SPtr<IntegrateValuesHelper> integrateValues;
    //!< compares velocity in integrateValues with target velocity and adjusts forcing accordingly.
    void collectData(double step);
    std::shared_ptr<vf::mpi::Communicator> comm;

private:
    double vx1Targed; //!< target velocity.
    double forcing;   //!< forcing at previous update step.
    double cellsVolume;
    double vx1Average;
    bool root;
    double Kpcrit; // Kp critical
    double Tcrit;  // the oscillation period
    double Tn;
    double Tv;
    double e;
    double Ta;
    double Kp;
    double Ki;
    double Kd;
    double y;
    double esum;
    double eold;
    // std::vector<CalcNodes> cnodes;
    std::string path;
};

#endif /* D3Q27RHODIFFERENCECoProcessor_H_ */
