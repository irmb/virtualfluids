/*
 *  TimeseriesPostprocessor.h
 *
 *  Created on: 08.05.2013
 *  Author: uphoff
 */

#ifndef TimeseriesPOSTPROCESSOR_H
#define TimeseriesPOSTPROCESSOR_H

#include "CoProcessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>

class TimeseriesCoProcessor;
typedef boost::shared_ptr<TimeseriesCoProcessor> TimeseriesCoProcessorPtr;

//! \brief     Writes timeseries of density and velocity to a file.
//! \details   Uses Integrate values helper, scheduler must be set in testcase.
//! \author    Sonja Uphoff
//! \date      May 2013

class TimeseriesCoProcessor: public CoProcessor {
public:
	TimeseriesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                                   D3Q27IntegrateValuesHelperPtr h1,
                                   const std::string& path, CommunicatorPtr comm);
	virtual ~TimeseriesCoProcessor();
	//! calls collectData.
   void process(double step); 
protected:
   //! object that can compute spacial average values in 3D-subdomain.
	D3Q27IntegrateValuesHelperPtr h1;  
	void collectData(double step);  
   CommunicatorPtr comm;
private:
	std::string path; //! output filename, e.g.  pathname + "/steps/timeseries"
   std::string fname;
};


#endif /* TimeseriesPOSTPROCESSOR_H */
