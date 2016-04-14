#include "Block3DConnectorFactory.h"
#include "D3Q27ETFullDirectConnector2.h"
#include "D3Q27ETFullVectorConnector.h"
#include "CoarseToFineNodeSetBlock3DConnector.h"
#include "FineToCoarseNodeSetBlock3DConnector.h"

Block3DConnectorFactory::Block3DConnectorFactory()
{
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorFactory::~Block3DConnectorFactory()
{
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorPtr Block3DConnectorFactory::createSameLevelDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir)
{
   return Block3DConnectorPtr(new D3Q27ETFullDirectConnector2(from, to, sendDir)); 
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorPtr Block3DConnectorFactory::createSameLevelVectorConnector(Block3DPtr block,
   VectorTransmitterPtr sender,
   VectorTransmitterPtr receiver,
   int sendDir)
{
   return Block3DConnectorPtr(new D3Q27ETFullVectorConnector(block, sender, receiver, sendDir));
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorPtr Block3DConnectorFactory::createCoarseToFineConnector(Block3DPtr block,
   VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
   VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
   VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
   VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
   int sendDir, D3Q27InterpolationProcessorPtr iprocessor)
{
   return Block3DConnectorPtr (new CoarseToFineNodeSetBlock3DConnector(block,
      sender00, receiver00, sender01, receiver01,
      sender10, receiver10, sender11, receiver11,
      sendDir, iprocessor));
}
//////////////////////////////////////////////////////////////////////////
Block3DConnectorPtr Block3DConnectorFactory::createFineToCoarseConnector(Block3DPtr block,
   VectorTransmitterPtr sender,
   VectorTransmitterPtr receiver,
   int sendDir,
   D3Q27InterpolationProcessorPtr iprocessor,
   FineToCoarseBlock3DConnector::CFconnectorType connType)
{
   return  Block3DConnectorPtr(new FineToCoarseNodeSetBlock3DConnector(block,
      sender, receiver, sendDir, iprocessor, connType));
}
