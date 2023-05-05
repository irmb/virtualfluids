#ifndef TimeDependentBCSimulationObserver_H
#define TimeDependentBCSimulationObserver_H

#include <PointerDefinitions.h>
#include <vector>

#include "SimulationObserver.h"

class Interactor3D;
class Grid3D;

//! \brief The class update interactors depend of time step.
//! \details TimeDependentBCSimulationObserver update every time step information in BCs throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCSimulationObserver : public SimulationObserver
{
public:
    TimeDependentBCSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
    ~TimeDependentBCSimulationObserver() override;

    void update(real step) override;

    //! add interactors to SimulationObserver
    void addInteractor(SPtr<Interactor3D> interactor);

private:
    std::vector<SPtr<Interactor3D>> interactors;
};

#endif /* TimeDependentBCSimulationObserver_H */
