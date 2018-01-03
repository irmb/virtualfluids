#include "CoProcessor.h"

#include <basics/utilities/UbScheduler.h>
#include <boost/bind.hpp>


CoProcessor::CoProcessor()
{
}

CoProcessor::CoProcessor(std::shared_ptr<Grid3D> grid, UbSchedulerPtr s): grid(grid), scheduler(s)
{
    connection = grid->connect(boost::bind(&CoProcessor::process, this, _1));
}

CoProcessor::~CoProcessor()
{
    grid->disconnect(connection);
}

void CoProcessor::disconnect()
{
    grid->disconnect(connection);
}

void CoProcessor::reconnect(std::shared_ptr<Grid3D> grid)
{
    this->grid = grid;
    this->grid->disconnect(connection);
    connection = this->grid->connect(boost::bind(&CoProcessor::process, this, _1));
}
