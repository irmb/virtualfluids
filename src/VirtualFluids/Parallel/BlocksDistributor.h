#ifndef BlocksDistributor_H
#define BlocksDistributor_H

#include "Communicator.h"
#include "Grid3D.h"

#include <memory>

class BlocksDistributor;
typedef std::shared_ptr<BlocksDistributor> BlocksDistributorPtr;

class BlocksDistributor
{
public:
   BlocksDistributor(Grid3DPtr grid, CommunicatorPtr comm);
   ~BlocksDistributor();

protected:
private:
   Grid3DPtr grid;
   CommunicatorPtr comm;
};

#endif


