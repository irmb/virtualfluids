#include "TimeDependentBCSimulationObserver.h"

#include "Grid3D.h"
#include "Interactor3D.h"
#include "UbScheduler.h"

TimeDependentBCSimulationObserver::TimeDependentBCSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s) : SimulationObserver(grid, s) {}
//////////////////////////////////////////////////////////////////////////
TimeDependentBCSimulationObserver::~TimeDependentBCSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCSimulationObserver::process(real step)
{
    if (scheduler->isDue(step)) {
        for (SPtr<Interactor3D> inter : interactors)
            inter->updateInteractor(step);
        UBLOG(logDEBUG3, "TimeDependentBCSimulationObserver::update:" << step);
    }
}
//////////////////////////////////////////////////////////////////////////
void TimeDependentBCSimulationObserver::addInteractor(SPtr<Interactor3D> interactor) { interactors.push_back(interactor); }

//////////////////////////////////////////////////////////////////////////
