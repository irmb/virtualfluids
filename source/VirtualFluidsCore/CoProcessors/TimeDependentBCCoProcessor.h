#ifndef TimeDependentBCPOSTPROCESSOR_H
#define TimeDependentBCPOSTPROCESSOR_H

#include "CoProcessor.h"
#include "Interactor3D.h"

#include <boost/shared_ptr.hpp>
class TimeDependentBCCoProcessor;
typedef boost::shared_ptr<TimeDependentBCCoProcessor> TimeDependentBCCoProcessorPtr;

//! \brief The class update interactors depend of time step. 
//! \details TimeDependentBCCoProcessor update every time step information in BCAdapters throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCCoProcessor: public CoProcessor {
public:
	TimeDependentBCCoProcessor(Grid3DPtr grid);
	virtual ~TimeDependentBCCoProcessor();
	void process(double step);
   //! add interactors to Postprocessor
   void addInteractor(Interactor3DPtr interactor);
protected:
private:
   std::vector<Interactor3DPtr> interactors;
};


#endif /* TimeDependentBCPOSTPROCESSOR_H */
