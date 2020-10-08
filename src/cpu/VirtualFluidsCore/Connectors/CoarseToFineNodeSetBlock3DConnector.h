//! \file CoarseToFineNodeSetBlock3DConnector.h
//! \class CoarseToFineNodeSetBlock3DConnector
//! \brief Connector interpolates and sends data from coarse level to fine.  
//! \author Konstantin Kutscher
//! \date 18.05.2015

#ifndef CoarseToFineNodeSetBlock3DConnector_H
#define CoarseToFineNodeSetBlock3DConnector_H

#include <vector>
#include <set>

//#include "basics/transmitter/TbTransmitter.h"
//#include "basics/transmitter/TbTransmitterLocal.h"
//#include "basics/container/CbVector.h"
#include "CoarseToFineBlock3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "InterpolationProcessor.h"
#include "MathUtil.hpp"
#include "Grid3D.h"
#include <PointerDefinitions.h>


class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)


// send direction:    E<->W     N<->S    T<->B  
//  ---------          x3       x3        x2    
// | 01 | 11 |         ^        ^         ^     
// |----+----|         +-> x2   +->x1     +->x1 
// | 00 | 10 |                                  
//  ---------                                   


class CoarseToFineNodeSetBlock3DConnector : public CoarseToFineBlock3DConnector
{
public:
   CoarseToFineNodeSetBlock3DConnector(SPtr<Block3D> block,
      VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
      VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
      VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
      VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
      int sendDir, InterpolationProcessorPtr iprocessor);

   void init() override;

   void fillSendVectors() override;
   void distributeReceiveVectors() override;

protected:
   typedef std::vector< int > INodeVector;
   using INodeSet = std::vector<INodeVector>;
   INodeSet  iNodeSetSender00;
   INodeSet  iNodeSetSender01;
   INodeSet  iNodeSetSender10;
   INodeSet  iNodeSetSender11;
   INodeSet  iNodeSetReceiver00;
   INodeSet  iNodeSetReceiver01;
   INodeSet  iNodeSetReceiver10;
   INodeSet  iNodeSetReceiver11;

   void writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF);
   void writeNodeToVector(vector_type& data, int& index, LBMReal* inode);
   void readICellCfromData(vector_type& data, int& index, LBMReal* icellC);

   void findCFCells();
   void findCFCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes);

   void findFCCells();
   void findFCCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes);

   void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);

   int bMaxX1, bMaxX2, bMaxX3;

   int minX1;
   int minX2;
   int minX3;

   int maxX1;
   int maxX2;
   int maxX3;

   int minHalfX1;
   int minHalfX2;
   int minHalfX3;

   int maxHalfX1;
   int maxHalfX2;
   int maxHalfX3;
};



#endif 
