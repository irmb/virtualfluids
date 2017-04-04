#ifndef D3Q27ADJUSTFORCINGCoProcessor_H
#define D3Q27ADJUSTFORCINGCoProcessor_H

#include "CoProcessor.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class AdjustForcingCoProcessor;
typedef boost::shared_ptr<AdjustForcingCoProcessor> AdjustForcingCoProcessorPtr;

//! \brief   Computes forcing such that a given velocity (vx1Targed) is reached inside an averaging domain (h1). 
//! \details Algorithm based on PID controller (proportional–integral–derivative controller). The parameters of PID controller estimation based on Ziegler–Nichols method. 
//!          Integrate values helper, scheduler must be set in test case.
//! \author: Konstantin Kutscher

class AdjustForcingCoProcessor: public CoProcessor {
public:
	AdjustForcingCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
                                   const std::string& path,
                                   IntegrateValuesHelperPtr integrateValues,
                                   double vTarged, CommunicatorPtr comm);
	virtual ~AdjustForcingCoProcessor();
	 //!< calls collect PostprocessData
   void process(double step);
protected:
   //!< object that can compute spacial average values in 3D-subdomain.
   IntegrateValuesHelperPtr integrateValues;
   //!< compares velocity in integrateValues with target velocity and adjusts forcing accordingly.
	void collectData(double step);  
   CommunicatorPtr comm;
private:
   double vx1Targed; //!< target velocity.
   double forcing; //!< forcing at previous update step. 
   double cellsVolume;
   double vx1Average;
   bool root;
   double Kpcrit; //Kp critical
   double Tcrit;  //the oscillation period 
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
   //std::vector<CalcNodes> cnodes;
   std::string path;
};


#endif /* D3Q27RHODIFFERENCECoProcessor_H_ */
