#ifndef NUMERICAL_TESTS_GRID_READER_H
#define NUMERICAL_TESTS_GRID_READER_H

#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"

#include <memory>

class InitialCondition;

class NumericalTestGridReader : public GridReader
{
public:
	static std::shared_ptr<NumericalTestGridReader> getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition);

protected:
	void setInitalNodeValues(const int numberOfNodes, const int level) const;
		
private:
	NumericalTestGridReader(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition);

	std::shared_ptr<InitialCondition> initialCondition;
	
};
#endif
