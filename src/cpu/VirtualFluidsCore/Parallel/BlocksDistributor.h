#ifndef BlocksDistributor_H
#define BlocksDistributor_H

#include <mpi/Communicator.h>
#include "Grid3D.h"

#include <PointerDefinitions.h>

class BlocksDistributor
{
public:
    BlocksDistributor(SPtr<Grid3D> grid, std::shared_ptr<vf::mpi::Communicator> comm);
    ~BlocksDistributor();

protected:
private:
    SPtr<Grid3D> grid;
    std::shared_ptr<vf::mpi::Communicator> comm;
};

#endif
