#ifndef ConnectorFactory_h__
#define ConnectorFactory_h__

#include "Block3DConnector.h"
#include "TransmitterType.h"
#include "InterpolationProcessor.h"
#include "FineToCoarseBlock3DConnector.h"

#include <PointerDefinitions.h>

class ConnectorFactory
{
public:
   ConnectorFactory() = default;;
   virtual ~ConnectorFactory() = default;;

   virtual SPtr<Block3DConnector> createSameLevelDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir) = 0;
   virtual SPtr<Block3DConnector> createSameLevelVectorConnector(SPtr<Block3D> block,
                                                              VectorTransmitterPtr sender, 
                                                              VectorTransmitterPtr receiver, 
                                                              int sendDir) = 0;
   virtual SPtr<Block3DConnector> createCoarseToFineConnector(SPtr<Block3D> block,
                                                            VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
                                                            VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
                                                            VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
                                                            VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
                                                            int sendDir, InterpolationProcessorPtr iprocessor) = 0;
   virtual SPtr<Block3DConnector> createFineToCoarseConnector(SPtr<Block3D> block, 
                                                           VectorTransmitterPtr sender, 
                                                           VectorTransmitterPtr receiver, 
                                                           int sendDir, 
                                                           InterpolationProcessorPtr iprocessor, 
                                                           FineToCoarseBlock3DConnector::CFconnectorType connType) = 0;

protected:
private:
};
#endif // ConnectorFactory_h__
