/*
 *  TimeseriesPostprocessor.h
 *
 *  Created on: 08.05.2013
 *  Author: uphoff
 */

#ifndef TimeseriesPOSTPROCESSOR_H
#define TimeseriesPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>

class TimeseriesPostprocessor;
typedef boost::shared_ptr<TimeseriesPostprocessor> TimeseriesWriterPostprocessorPtr;

//! \brief     Writes timeseries of density and velocity to a file.
//! \details   Uses Integrate values helper, scheduler must be set in testcase.
//! \author    Sonja Uphoff
//! \date      May 2013

class TimeseriesPostprocessor: public Postprocessor {
public:
	TimeseriesPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                                   D3Q27IntegrateValuesHelperPtr h1,
                                   const std::string& path, CommunicatorPtr comm);
	virtual ~TimeseriesPostprocessor();
	//! calls collect PostprocessData.
   void update(double step); 
protected:
   //! object that can compute spacial average values in 3D-subdomain.
	D3Q27IntegrateValuesHelperPtr h1;  
	void collectPostprocessData(double step);  
   CommunicatorPtr comm;
private:
	std::string path; //! output filename, e.g.  pathname + "/steps/timeseries"
   std::string fname;
};


#endif /* TimeseriesPOSTPROCESSOR_H */
