/*
 *  TimeseriesSimulationObserver.h
 *
 *  Created on: 08.05.2013
 *  Author: uphoff
 */

#ifndef TimeseriesSimulationObserver_H
#define TimeseriesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class IntegrateValuesHelper;

//! \brief     Writes timeseries of density and velocity to a file.
//! \details   Uses Integrate values helper, scheduler must be set in testcase.
//! \author    Sonja Uphoff
//! \date      May 2013

class TimeseriesSimulationObserver : public SimulationObserver
{
public:
    TimeseriesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<IntegrateValuesHelper> h1,
                          const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm);
    ~TimeseriesSimulationObserver() override;

    //! calls collectData.
    void update(real step) override;

protected:
    void collectData(real step);

    //! object that can compute spacial average values in 3D-subdomain.
    SPtr<IntegrateValuesHelper> h1;
    std::shared_ptr<vf::mpi::Communicator> comm;

private:
    std::string path; //! output filename, e.g.  pathname + "/steps/timeseries"
    std::string fname;
};

#endif /* TimeseriesSimulationObserver_H */
