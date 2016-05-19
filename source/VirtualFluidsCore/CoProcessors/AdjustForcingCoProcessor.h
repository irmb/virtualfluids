#ifndef D3Q27ADJUSTFORCINGCoProcessor_H
#define D3Q27ADJUSTFORCINGCoProcessor_H

#include "CoProcessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class AdjustForcingCoProcessor;
typedef boost::shared_ptr<AdjustForcingCoProcessor> AdjustForcingCoProcessorPtr;

//! \brief   Computes Forcing such that a given velocity (vxZiel) is reached inside an averaging domain (h1). 
//! \details Integrate values helper, scheduler must be set in test case. Example usage: bKanal.cpp
//! \author: Konstantin Kutscher

class AdjustForcingCoProcessor: public CoProcessor {
public:
	AdjustForcingCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                                   const std::string& path,
                                   D3Q27IntegrateValuesHelperPtr integrateValues,
                                   double vTarged, double forcing, CommunicatorPtr comm);
	virtual ~AdjustForcingCoProcessor();
	 //!< calls collect PostprocessData
   void process(double step);
protected:
   //!< object that can compute spacial average values in 3D-subdomain.
   D3Q27IntegrateValuesHelperPtr integrateValues;
   //!< compares velocity in integrateValues with target velocity and adjusts forcing accordingly.
	void collectData(double step);  
   CommunicatorPtr comm;
private:
   double vTarged; //!< target velocity.
   double forcing; //!< forcing at previous update step. 
   std::vector<CalcNodes> cnodes;
   std::string path;
};


#endif /* D3Q27RHODIFFERENCECoProcessor_H_ */
