#ifndef TimeDependentBCCoProcessor_H
#define TimeDependentBCCoProcessor_H

#include <PointerDefinitions.h>
#include <vector>

#include "CoProcessor.h"

class Interactor3D;
class Grid3D;

//! \brief The class update interactors depend of time step.
//! \details TimeDependentBCCoProcessor update every time step information in BCs throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCCoProcessor : public CoProcessor
{
public:
    TimeDependentBCCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
    ~TimeDependentBCCoProcessor() override;

    void process(real step) override;

    //! add interactors to CoProcessor
    void addInteractor(SPtr<Interactor3D> interactor);

private:
    std::vector<SPtr<Interactor3D>> interactors;
};

#endif /* TimeDependentBCCoProcessor_H */
