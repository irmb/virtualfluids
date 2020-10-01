#ifndef Block3DConnectorFactory_h__
#define Block3DConnectorFactory_h__

#include "ConnectorFactory.h"

#include <PointerDefinitions.h>

class Block3DConnectorFactory : public ConnectorFactory
{
public:
   Block3DConnectorFactory();
   virtual ~Block3DConnectorFactory();

   virtual SPtr<Block3DConnector> createSameLevelDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);

   virtual SPtr<Block3DConnector> createSameLevelVectorConnector(SPtr<Block3D> block,
      VectorTransmitterPtr sender,
      VectorTransmitterPtr receiver,
      int sendDir);

   virtual SPtr<Block3DConnector> createCoarseToFineConnector(SPtr<Block3D> block,
      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
      int sendDir, InterpolationProcessorPtr iprocessor);

   virtual SPtr<Block3DConnector> createFineToCoarseConnector(SPtr<Block3D> block,
      VectorTransmitterPtr sender,
      VectorTransmitterPtr receiver,
      int sendDir,
      InterpolationProcessorPtr iprocessor,
      FineToCoarseBlock3DConnector::CFconnectorType connType);

private:

};
#endif // Block3DConnectorFactory_h__

