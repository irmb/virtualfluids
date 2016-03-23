#ifndef D3Q27ADJUSTFORCINGPOSTPROCESSOR_H
#define D3Q27ADJUSTFORCINGPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"

#include <boost/shared_ptr.hpp>
class D3Q27AdjustForcingPostprocessor;
typedef boost::shared_ptr<D3Q27AdjustForcingPostprocessor> D3Q27AdjustForcingPostprocessorPtr;

//! \brief   Computes Forcing such that a given velocity (vxZiel) is reached inside an averaging domain (h1). 
//! \details Integrate values helper, scheduler must be set in test case. Example usage: bKanal.cpp
//! \author: Konstantin Kutscher

class D3Q27AdjustForcingPostprocessor: public Postprocessor {
public:
	D3Q27AdjustForcingPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                                   const std::string& path,
                                   D3Q27IntegrateValuesHelperPtr integrateValues,
                                   LBMReal vTarged, LBMReal forcing, CommunicatorPtr comm);
	virtual ~D3Q27AdjustForcingPostprocessor();
	 //!< calls collect PostprocessData
   void update(double step);
protected:
   //!< object that can compute spacial average values in 3D-subdomain.
   D3Q27IntegrateValuesHelperPtr integrateValues;
   //!< compares velocity in integrateValues with target velocity and adjusts forcing accordingly.
	void collectPostprocessData(double step);  
   CommunicatorPtr comm;
private:
   double vPreviousStep; //!< velocity at previous update step.
   double vTarged; //!< target velocity.
   double forcing; //!< forcing at previous update step. 
   std::vector<CalcNodes> cnodes;
   std::string path;
};


#endif /* D3Q27RHODIFFERENCEPOSTPROCESSOR_H_ */
