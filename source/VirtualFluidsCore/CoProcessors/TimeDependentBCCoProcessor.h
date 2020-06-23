#ifndef TimeDependentBCCoProcessor_H
#define TimeDependentBCCoProcessor_H

#include <vector>
#include <PointerDefinitions.h>

#include "CoProcessor.h"

class Interactor3D;
class Grid3D;

//! \brief The class update interactors depend of time step. 
//! \details TimeDependentBCCoProcessor update every time step information in BCAdapters throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCCoProcessor : public CoProcessor
{
public:
	TimeDependentBCCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s);
	virtual ~TimeDependentBCCoProcessor();

	void process(double step) override;

   //! add interactors to CoProcessor
   void addInteractor(SPtr<Interactor3D> interactor);

private:
   std::vector<SPtr<Interactor3D> > interactors;
};


#endif /* TimeDependentBCCoProcessor_H */
