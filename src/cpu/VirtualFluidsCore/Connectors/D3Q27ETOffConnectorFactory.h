//#ifndef D3Q27ETOffConnectorFactory_h__
//#define D3Q27ETOffConnectorFactory_h__
//
//#include "Block3DConnectorFactory.h"
//
//#include <PointerDefinitions.h>
// class D3Q27ETOffConnectorFactory;
// typedef SPtr<D3Q27ETOffConnectorFactory> D3Q27ETOffSPtr<ConnectorFactory>;
//
// class D3Q27ETOffConnectorFactory : public Block3DConnectorFactory
//{
// public:
//   D3Q27ETOffConnectorFactory();
//   virtual ~D3Q27ETOffConnectorFactory();
//
//   virtual SPtr<Block3DConnector> createCoarseToFineConnector(SPtr<Block3D> block,
//      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
//      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
//      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
//      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
//      int sendDir, D3Q27InterpolationProcessorPtr iprocessor);
//
//   virtual SPtr<Block3DConnector> createFineToCoarseConnector(SPtr<Block3D> block,
//      VectorTransmitterPtr sender,
//      VectorTransmitterPtr receiver,
//      int sendDir,
//      D3Q27InterpolationProcessorPtr iprocessor,
//      FineToCoarseBlock3DConnector::CFconnectorType connType);
//
// private:
//
//};
//#endif // D3Q27ETOffConnectorFactory_h__
