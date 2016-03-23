#ifndef D3Q27ETFULLVECTORCONNECTOR_H
#define D3Q27ETFULLVECTORCONNECTOR_H

#include <vector>

#include "RemoteBlock3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernelETD3Q27.h"
#include "EsoTwistD3Q27System.h"


//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)
class D3Q27ETFullVectorConnector : public RemoteBlock3DConnector
{
public:
   D3Q27ETFullVectorConnector(  Block3DPtr block
      , VectorTransmitterPtr sender
      , VectorTransmitterPtr receiver
      , int sendDir); 

   void init();

   void fillSendVectors();
   void distributeReceiveVectors();

protected:
   void fillData(EsoTwist3DPtr  fFrom, int& index, const int& x1, const int& x2, const int& x3);
   void distributeData(EsoTwist3DPtr  fTo, int& index, const int& x1, const int& x2, const int& x3);

};



#endif //D3Q27VECTORCONNECTOR_H

