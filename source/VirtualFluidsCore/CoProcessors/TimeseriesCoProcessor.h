/*
 *  TimeseriesCoProcessor.h
 *
 *  Created on: 08.05.2013
 *  Author: uphoff
 */

#ifndef TimeseriesCoProcessor_H
#define TimeseriesCoProcessor_H

#include <memory>
#include <string>

#include "CoProcessor.h"

class Communicator;
class Grid3D;
class UbScheduler;
class IntegrateValuesHelper;

class TimeseriesCoProcessor;
typedef std::shared_ptr<TimeseriesCoProcessor> TimeseriesCoProcessorPtr;

//! \brief     Writes timeseries of density and velocity to a file.
//! \details   Uses Integrate values helper, scheduler must be set in testcase.
//! \author    Sonja Uphoff
//! \date      May 2013

class TimeseriesCoProcessor : public CoProcessor
{
public:
    TimeseriesCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<IntegrateValuesHelper> h1, const std::string& path, std::shared_ptr<Communicator> comm);
    virtual ~TimeseriesCoProcessor();

    //! calls collectData.
    void process(double step) override;

protected:
    void collectData(double step);

    //! object that can compute spacial average values in 3D-subdomain.
    std::shared_ptr<IntegrateValuesHelper> h1;
    std::shared_ptr<Communicator> comm;

private:
    std::string path; //! output filename, e.g.  pathname + "/steps/timeseries"
    std::string fname;
};


#endif /* TimeseriesCoProcessor_H */
