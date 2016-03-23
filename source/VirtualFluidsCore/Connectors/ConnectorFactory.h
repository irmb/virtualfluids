#ifndef ConnectorFactory_h__
#define ConnectorFactory_h__

#include "Block3DConnector.h"
#include "TransmitterType.h"
#include "D3Q27InterpolationProcessor.h"
#include "FineToCoarseBlock3DConnector.h"

#include <boost/shared_ptr.hpp>
class ConnectorFactory;
typedef boost::shared_ptr<ConnectorFactory> ConnectorFactoryPtr;

class ConnectorFactory
{
public:
   ConnectorFactory() {};
   virtual ~ConnectorFactory() {};

   virtual Block3DConnectorPtr createSameLevelDirectConnector(Block3DPtr from, Block3DPtr to, int sendDir) = 0;
   virtual Block3DConnectorPtr createSameLevelVectorConnector(Block3DPtr block,
                                                              VectorTransmitterPtr sender, 
                                                              VectorTransmitterPtr receiver, 
                                                              int sendDir) = 0;
   virtual Block3DConnectorPtr createCoarseToFineConnector(Block3DPtr block,
                                                            VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
                                                            VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
                                                            VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
                                                            VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
                                                            int sendDir, D3Q27InterpolationProcessorPtr iprocessor) = 0;
   virtual Block3DConnectorPtr createFineToCoarseConnector(Block3DPtr block, 
                                                           VectorTransmitterPtr sender, 
                                                           VectorTransmitterPtr receiver, 
                                                           int sendDir, 
                                                           D3Q27InterpolationProcessorPtr iprocessor, 
                                                           FineToCoarseBlock3DConnector::CFconnectorType connType) = 0;

protected:
private:
};
#endif // ConnectorFactory_h__
