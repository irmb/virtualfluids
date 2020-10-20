//#include "D3Q27ETOffConnectorFactory.h"
//#include "TransmitterType.h"
//#include "D3Q27ETCFOffVectorConnector.h"
//#include "D3Q27ETFCOffVectorConnector.h"
//#include "D3Q27ETFCVectorConnector.h"
//#include "FineToCoarseBlock3DConnector.h"
//
// D3Q27ETOffConnectorFactory::D3Q27ETOffConnectorFactory()
//{
//}
////////////////////////////////////////////////////////////////////////////
// D3Q27ETOffConnectorFactory::~D3Q27ETOffConnectorFactory()
//{
//}
////////////////////////////////////////////////////////////////////////////
// SPtr<Block3DConnector> D3Q27ETOffConnectorFactory::createCoarseToFineConnector(SPtr<Block3D> block,
//   VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
//   VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
//   VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
//   VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
//   int sendDir, D3Q27InterpolationProcessorPtr iprocessor)
//{
//   return SPtr<Block3DConnector>(new D3Q27ETCFOffVectorConnector<VectorTransmitter>(block,
//      sender00, receiver00, sender01, receiver01,
//      sender10, receiver10, sender11, receiver11,
//      sendDir, iprocessor));
//}
////////////////////////////////////////////////////////////////////////////
// SPtr<Block3DConnector> D3Q27ETOffConnectorFactory::createFineToCoarseConnector(SPtr<Block3D> block,
//   VectorTransmitterPtr sender,
//   VectorTransmitterPtr receiver,
//   int sendDir,
//   D3Q27InterpolationProcessorPtr iprocessor,
//   FineToCoarseBlock3DConnector::CFconnectorType connType)
//{
//   return  SPtr<Block3DConnector>(new D3Q27ETFCOffVectorConnector<VectorTransmitter>(block,
//      sender, receiver, sendDir, iprocessor, connType));
//}
