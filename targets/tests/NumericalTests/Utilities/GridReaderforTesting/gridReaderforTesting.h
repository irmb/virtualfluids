#ifndef GRIDREADER_FOR_TESTING_H
#define GRIDREADER_FOR_TESTING_H

#include "VirtualFluids_GPU/DataStructureInitializer/GridReaderFiles/GridReader.h"

#include <memory>

class InitialCondition;

class GridReaderforTesting : public GridReader 
{
public:
	
	void setInitalNodeValues(const int numberOfNodes, const int level) const;

	GridReaderforTesting(std::shared_ptr<Parameter> para, std::shared_ptr<InitialCondition> initialCondition);
			
private:
	std::shared_ptr<InitialCondition> initialCondition;
	
};
#endif
