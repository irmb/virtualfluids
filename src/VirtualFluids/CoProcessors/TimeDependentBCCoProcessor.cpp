#include "TimeDependentBCCoProcessor.h"
#include <boost/foreach.hpp>

using namespace std;

TimeDependentBCCoProcessor::TimeDependentBCCoProcessor(Grid3DPtr grid) : CoProcessor(grid,  UbSchedulerPtr(new UbScheduler(1)))
{

}
//////////////////////////////////////////////////////////////////////////
TimeDependentBCCoProcessor::~TimeDependentBCCoProcessor() 
{
	
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::process(double step)
{
   BOOST_FOREACH(Interactor3DPtr inter, interactors)
      inter->updateInteractor( step );
   UBLOG(logDEBUG3, "TimeDependentBCCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::addInteractor( Interactor3DPtr interactor )
{
   interactors.push_back(interactor);
}

//////////////////////////////////////////////////////////////////////////


