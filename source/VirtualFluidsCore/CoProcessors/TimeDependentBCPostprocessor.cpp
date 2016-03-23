#include "TimeDependentBCPostprocessor.h"
#include <boost/foreach.hpp>

using namespace std;

TimeDependentBCPostprocessor::TimeDependentBCPostprocessor(Grid3DPtr grid) : Postprocessor(grid,  UbSchedulerPtr(new UbScheduler(1)))
{

}
//////////////////////////////////////////////////////////////////////////
TimeDependentBCPostprocessor::~TimeDependentBCPostprocessor() 
{
	
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCPostprocessor::update(double step)
{
   BOOST_FOREACH(Interactor3DPtr inter, interactors)
      inter->updateInteractor( step );
   UBLOG(logDEBUG3, "TimeDependentBCPostprocessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCPostprocessor::addInteractor( Interactor3DPtr interactor )
{
   interactors.push_back(interactor);
}

//////////////////////////////////////////////////////////////////////////


