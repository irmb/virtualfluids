#include "TimeDependentBCCoProcessor.h"

#include "Grid3D.h"
#include "Interactor3D.h"
#include "UbScheduler.h"

TimeDependentBCCoProcessor::TimeDependentBCCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s) : CoProcessor(grid, s) {}
//////////////////////////////////////////////////////////////////////////
TimeDependentBCCoProcessor::~TimeDependentBCCoProcessor() = default;
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::process(real step)
{
    if (scheduler->isDue(step)) {
        for (SPtr<Interactor3D> inter : interactors)
            inter->updateInteractor(step);
        UBLOG(logDEBUG3, "TimeDependentBCCoProcessor::update:" << step);
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCCoProcessor::addInteractor(SPtr<Interactor3D> interactor) { interactors.push_back(interactor); }

//////////////////////////////////////////////////////////////////////////
