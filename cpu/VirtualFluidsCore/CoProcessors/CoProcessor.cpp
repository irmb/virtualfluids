#include "CoProcessor.h"

#include "Grid3D.h"
#include "UbScheduler.h"

CoProcessor::CoProcessor()
{
}

CoProcessor::CoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s): grid(grid), scheduler(s)
{

}

CoProcessor::~CoProcessor()
{

}

