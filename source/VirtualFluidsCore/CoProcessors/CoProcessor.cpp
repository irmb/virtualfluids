#include "CoProcessor.h"

#include "Grid3D.h"
#include "UbScheduler.h"

CoProcessor::CoProcessor()
{
}

CoProcessor::CoProcessor(std::shared_ptr<Grid3D> grid, UbSchedulerPtr s): grid(grid), scheduler(s)
{

}

CoProcessor::~CoProcessor()
{

}

