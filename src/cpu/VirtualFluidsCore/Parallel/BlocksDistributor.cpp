#include "BlocksDistributor.h"

BlocksDistributor::BlocksDistributor(SPtr<Grid3D> grid, SPtr<Communicator> comm) : grid(grid), comm(comm) {}
