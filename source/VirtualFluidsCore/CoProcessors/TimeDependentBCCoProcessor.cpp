#include "TimeDependentBCCoProcessor.h"

#include "Interactor3D.h"
#include "UbScheduler.h"
#include "Grid3D.h"

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
   for(Interactor3DPtr inter : interactors)
      inter->updateInteractor( step );
   UBLOG(logDEBUG3, "TimeDependentBCCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::addInteractor( Interactor3DPtr interactor )
{
   interactors.push_back(interactor);
}

//////////////////////////////////////////////////////////////////////////


