#ifndef GRIDREADER_FOR_TESTING_H
#define GRIDREADER_FOR_TESTING_H

#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"

#include <memory>

class InitialCondition;

class GridReaderforTesting : public GridReader 
{
public:
	static std::shared_ptr<GridReaderforTesting> getNewInstance(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition);
	void setInitalNodeValues(const int numberOfNodes, const int level) const;
		
private:
	GridReaderforTesting(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition);

	std::shared_ptr<InitialCondition> initialCondition;
	
};
#endif
