#ifndef TimeDependentBCPOSTPROCESSOR_H
#define TimeDependentBCPOSTPROCESSOR_H

#include "Postprocessor.h"
#include "Interactor3D.h"

#include <boost/shared_ptr.hpp>
class TimeDependentBCPostprocessor;
typedef boost::shared_ptr<TimeDependentBCPostprocessor> TimeDependentBCPostprocessorPtr;

//! \brief The class update interactors depend of time step. 
//! \details TimeDependentBCPostprocessor update every time step information in BCAdapters throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCPostprocessor: public Postprocessor {
public:
	TimeDependentBCPostprocessor(Grid3DPtr grid);
	virtual ~TimeDependentBCPostprocessor();
	void update(double step);
   //! add interactors to Postprocessor
   void addInteractor(Interactor3DPtr interactor);
protected:
private:
   std::vector<Interactor3DPtr> interactors;
};


#endif /* TimeDependentBCPOSTPROCESSOR_H */
