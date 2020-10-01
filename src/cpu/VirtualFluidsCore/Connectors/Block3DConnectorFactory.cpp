#include "Block3DConnectorFactory.h"
#include "D3Q27ETFullDirectConnector.h"
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
SPtr<Block3DConnector> Block3DConnectorFactory::createSameLevelDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir)
{
   return SPtr<Block3DConnector>(new D3Q27ETFullDirectConnector(from, to, sendDir)); 
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3DConnector> Block3DConnectorFactory::createSameLevelVectorConnector(SPtr<Block3D> block,
   VectorTransmitterPtr sender,
   VectorTransmitterPtr receiver,
   int sendDir)
{
   return SPtr<Block3DConnector>(new D3Q27ETFullVectorConnector(block, sender, receiver, sendDir));
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3DConnector> Block3DConnectorFactory::createCoarseToFineConnector(SPtr<Block3D> block,
   VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
   VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
   VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
   VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
   int sendDir, InterpolationProcessorPtr iprocessor)
{
   return SPtr<Block3DConnector> (new CoarseToFineNodeSetBlock3DConnector(block,
      sender00, receiver00, sender01, receiver01,
      sender10, receiver10, sender11, receiver11,
      sendDir, iprocessor));
}
//////////////////////////////////////////////////////////////////////////
SPtr<Block3DConnector> Block3DConnectorFactory::createFineToCoarseConnector(SPtr<Block3D> block,
   VectorTransmitterPtr sender,
   VectorTransmitterPtr receiver,
   int sendDir,
   InterpolationProcessorPtr iprocessor,
   FineToCoarseBlock3DConnector::CFconnectorType connType)
{
   return  SPtr<Block3DConnector>(new FineToCoarseNodeSetBlock3DConnector(block,
      sender, receiver, sendDir, iprocessor, connType));
}
