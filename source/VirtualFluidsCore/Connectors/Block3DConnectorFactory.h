#ifndef Block3DConnectorFactory_h__
#define Block3DConnectorFactory_h__

#include "ConnectorFactory.h"

#include <memory>
class Block3DConnectorFactory;
typedef std::shared_ptr<Block3DConnectorFactory> Block3DConnectorFactoryPtr;

class Block3DConnectorFactory : public ConnectorFactory
{
public:
   Block3DConnectorFactory();
   virtual ~Block3DConnectorFactory();

   virtual Block3DConnectorPtr createSameLevelDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir);

   virtual Block3DConnectorPtr createSameLevelVectorConnector(Block3DPtr block,
      VectorTransmitterPtr sender,
      VectorTransmitterPtr receiver,
      int sendDir);

   virtual Block3DConnectorPtr createCoarseToFineConnector(Block3DPtr block,
      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
      int sendDir, InterpolationProcessorPtr iprocessor);

   virtual Block3DConnectorPtr createFineToCoarseConnector(Block3DPtr block,
      VectorTransmitterPtr sender,
      VectorTransmitterPtr receiver,
      int sendDir,
      InterpolationProcessorPtr iprocessor,
      FineToCoarseBlock3DConnector::CFconnectorType connType);

private:

};
#endif // Block3DConnectorFactory_h__

