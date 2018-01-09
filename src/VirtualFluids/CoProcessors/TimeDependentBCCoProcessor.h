#ifndef TimeDependentBCCoProcessor_H
#define TimeDependentBCCoProcessor_H

#include <vector>
#include <memory>

#include "CoProcessor.h"

class Interactor3D;
class Grid3D;

class TimeDependentBCCoProcessor;
typedef std::shared_ptr<TimeDependentBCCoProcessor> TimeDependentBCCoProcessorPtr;

//! \brief The class update interactors depend of time step. 
//! \details TimeDependentBCCoProcessor update every time step information in BCAdapters throw Interactors
//! \author Sonja Uphoff, Kostyantyn Kucher
class TimeDependentBCCoProcessor : public CoProcessor
{
public:
	TimeDependentBCCoProcessor(std::shared_ptr<Grid3D> grid);
	virtual ~TimeDependentBCCoProcessor();

	void process(double step) override;

   //! add interactors to CoProcessor
   void addInteractor(std::shared_ptr<Interactor3D> interactor);

private:
   std::vector<std::shared_ptr<Interactor3D> > interactors;
};


#endif /* TimeDependentBCCoProcessor_H */
