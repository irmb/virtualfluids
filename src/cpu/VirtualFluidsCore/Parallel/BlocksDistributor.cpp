#include "BlocksDistributor.h"

BlocksDistributor::BlocksDistributor(SPtr<Grid3D> grid, std::shared_ptr<vf::mpi::Communicator> comm) : grid(grid), comm(comm) {}
