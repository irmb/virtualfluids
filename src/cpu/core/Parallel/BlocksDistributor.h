#ifndef BlocksDistributor_H
#define BlocksDistributor_H

#include <parallel/Communicator.h>
#include "Grid3D.h"

#include <PointerDefinitions.h>

class BlocksDistributor
{
public:
    BlocksDistributor(SPtr<Grid3D> grid, std::shared_ptr<vf::parallel::Communicator> comm);
    ~BlocksDistributor();

protected:
private:
    SPtr<Grid3D> grid;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif
