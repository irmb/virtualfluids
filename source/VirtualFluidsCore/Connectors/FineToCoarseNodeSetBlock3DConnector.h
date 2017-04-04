/// \file CoarseToFineNodeSetBlock3DConnector.h
/// \class CoarseToFineNodeSetBlock3DConnector
/// \brief Connector interpolates and sends data from coarse level to fine.  
/// \author Konstantin Kutscher
/// \date 21.05.2015

#ifndef FineToCoarseNodeSetBlock3DConnector_H
#define FineToCoarseNodeSetBlock3DConnector_H

#include <vector>
#include "FineToCoarseBlock3DConnector.h"
#include "D3Q27System.h"
#include "Block3D.h"
#include "Grid3D.h"
#include "LBMKernel.h"
#include "InterpolationProcessor.h"
#include "MathUtil.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>

class Block3D;

//daten werden in einen vector (dieser befindet sich im transmitter) kopiert
//der vector wird via transmitter uebertragen
//transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
//transmitter sein, der von Transmitter abgeleitet ist ;-)

class FineToCoarseNodeSetBlock3DConnector : public FineToCoarseBlock3DConnector
{
public:
   FineToCoarseNodeSetBlock3DConnector(Block3DPtr block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver, int sendDir, InterpolationProcessorPtr iprocessor, CFconnectorType connType);
   void init();
   void fillSendVectors();
   void distributeReceiveVectors();
protected:
   typedef std::vector< int > INodeVector;
   typedef std::vector < INodeVector > INodeSet;
   INodeSet  iNodeSetSender;
   INodeSet  iNodeSetReceiver;

   void readICellFfromData(vector_type& data, int& index, D3Q27ICell& icellF);
   void readNodeFromVector(vector_type& data, int& index, LBMReal* inode);

   void writeICellCtoData(vector_type& data, int& index, LBMReal* icellC);

   void findFCCells();
   void findFCCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes);

   void findCFCells();
   void findCFCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes);

   //void getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3);


   int bMaxX1, bMaxX2, bMaxX3;
   
   int minX1;
   int minX2;
   int minX3;

   int maxX1;
   int maxX2;
   int maxX3;

   int minOffX1;
   int minOffX2;
   int minOffX3;

   int maxOffX1;
   int maxOffX2;
   int maxOffX3;
};





#endif
