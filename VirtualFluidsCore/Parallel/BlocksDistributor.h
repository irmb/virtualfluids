#ifndef BlocksDistributor_H
#define BlocksDistributor_H

#include "Communicator.h"
#include "Grid3D.h"

#include <boost/shared_ptr.hpp>

class BlocksDistributor;
typedef boost::shared_ptr<BlocksDistributor> BlocksDistributorPtr;

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


