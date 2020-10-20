#ifndef TransmitterType_h__
#define TransmitterType_h__

#include "D3Q27System.h"
#include "basics/container/CbVector.h"
#include "basics/transmitter/TbTransmitter.h"
#include "basics/transmitter/TbTransmitterLocal.h"
#include <PointerDefinitions.h>

using VectorTransmitter    = TbTransmitter<CbVector<LBMReal>>;
using vector_type          = VectorTransmitter::value_type;
using VectorTransmitterPtr = SPtr<TbTransmitter<CbVector<LBMReal>>>;

#endif // TransmitterType_h__
