#ifndef BlocksDistributor_H
#define BlocksDistributor_H

#include "Communicator.h"
#include "Grid3D.h"

#include <PointerDefinitions.h>

class BlocksDistributor
{
public:
    BlocksDistributor(SPtr<Grid3D> grid, SPtr<Communicator> comm);
    ~BlocksDistributor();

protected:
private:
    SPtr<Grid3D> grid;
    SPtr<Communicator> comm;
};

#endif
