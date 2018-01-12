#include "TimeDependentBCCoProcessor.h"

#include "Interactor3D.h"
#include "UbScheduler.h"
#include "Grid3D.h"

TimeDependentBCCoProcessor::TimeDependentBCCoProcessor(SPtr<Grid3D> grid) : CoProcessor(grid,  SPtr<UbScheduler>(new UbScheduler(1)))
{

}
//////////////////////////////////////////////////////////////////////////
TimeDependentBCCoProcessor::~TimeDependentBCCoProcessor() 
{
	
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::process(double step)
{
   for(SPtr<Interactor3D> inter : interactors)
      inter->updateInteractor( step );
   UBLOG(logDEBUG3, "TimeDependentBCCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::addInteractor( SPtr<Interactor3D> interactor )
{
   interactors.push_back(interactor);
}

//////////////////////////////////////////////////////////////////////////


