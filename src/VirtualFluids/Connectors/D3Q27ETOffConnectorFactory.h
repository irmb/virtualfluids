//#ifndef D3Q27ETOffConnectorFactory_h__
//#define D3Q27ETOffConnectorFactory_h__
//
//#include "Block3DConnectorFactory.h"
//
//#include <memory>
//class D3Q27ETOffConnectorFactory;
//typedef std::shared_ptr<D3Q27ETOffConnectorFactory> D3Q27ETOffConnectorFactoryPtr;
//
//class D3Q27ETOffConnectorFactory : public Block3DConnectorFactory
//{
//public:
//   D3Q27ETOffConnectorFactory();
//   virtual ~D3Q27ETOffConnectorFactory();
//
//   virtual Block3DConnectorPtr createCoarseToFineConnector(Block3DPtr block,
//      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
//      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
//      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
//      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
//      int sendDir, D3Q27InterpolationProcessorPtr iprocessor);
//
//   virtual Block3DConnectorPtr createFineToCoarseConnector(Block3DPtr block,
//      VectorTransmitterPtr sender,
//      VectorTransmitterPtr receiver,
//      int sendDir,
//      D3Q27InterpolationProcessorPtr iprocessor,
//      FineToCoarseBlock3DConnector::CFconnectorType connType);
//
//private:
//
//};
//#endif // D3Q27ETOffConnectorFactory_h__

