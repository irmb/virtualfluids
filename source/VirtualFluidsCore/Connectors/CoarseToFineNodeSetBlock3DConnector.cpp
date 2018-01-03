#include "CoarseToFineNodeSetBlock3DConnector.h"
#include "DataSet3D.h"


////////////////////////////////////////////////////////////////////////////
CoarseToFineNodeSetBlock3DConnector::CoarseToFineNodeSetBlock3DConnector(Block3DPtr block,
   VectorTransmitterPtr sender00, VectorTransmitterPtr receiver00,
   VectorTransmitterPtr sender01, VectorTransmitterPtr receiver01,
   VectorTransmitterPtr sender10, VectorTransmitterPtr receiver10,
   VectorTransmitterPtr sender11, VectorTransmitterPtr receiver11,
   int sendDir, InterpolationProcessorPtr iprocessor) : CoarseToFineBlock3DConnector(block, sender00, receiver00,
   sender01, receiver01,
   sender10, receiver10,
   sender11, receiver11,
   sendDir, iprocessor)
{
}
//////////////////////////////////////////////////////////////////////////
void CoarseToFineNodeSetBlock3DConnector::init()
{
   bMaxX1 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1();
   bMaxX2 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2();
   bMaxX3 = (int)block.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3();

   minX1 = 0;
   minX2 = 0;
   minX3 = 0;
   maxX1 = bMaxX1 - 1;
   maxX2 = bMaxX2 - 1;
   maxX3 = bMaxX3 - 1;

   minHalfX1 = 0;
   minHalfX2 = 0;
   minHalfX3 = 0;

   maxHalfX1 = 0;
   maxHalfX2 = 0;
   maxHalfX3 = 0;

   if (Utilities::isEven(bMaxX1))
   {
      minHalfX1 = bMaxX1 / 2 - 1;
      maxHalfX1 = bMaxX1 / 2 - 1;
   }
   else if (Utilities::isOdd(bMaxX1))
   {
      minHalfX1 = bMaxX1 / 2;
      maxHalfX1 = bMaxX1 / 2 - 1;
   }

   if (Utilities::isEven(bMaxX2))
   {
      minHalfX2 = bMaxX2 / 2 - 1;
      maxHalfX2 = bMaxX2 / 2 - 1;
   }
   else if (Utilities::isOdd(bMaxX2))
   {
      minHalfX2 = bMaxX2 / 2;
      maxHalfX2 = bMaxX2 / 2 - 1;
   }

   if (Utilities::isEven(bMaxX3))
   {
      minHalfX3 = bMaxX3 / 2 - 1;
      maxHalfX3 = bMaxX3 / 2 - 1;
   }
   else if (Utilities::isOdd(bMaxX3))
   {
      minHalfX3 = bMaxX3 / 2;
      maxHalfX3 = bMaxX3 / 2 - 1;
   }

   //int       sendSize = 0;
   LBMReal initValue = -999.0;

   int sendDataPerNode = 27/*f*/;
   int iCellSize = 8; //size of interpolation cell

   findCFCells();
   findFCCells();

   //////////////////////////////////////////////////////
   //Debug
   //////////////////////////////////////////////////////
   if (block.lock()->getGlobalID() == 2234)
   {
      int test = 0;
   }

   if (sender00) sender00->getData().resize(iNodeSetSender00.size()*iCellSize*sendDataPerNode, initValue);
   else sender00 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (sender01)  sender01->getData().resize(iNodeSetSender01.size()*iCellSize*sendDataPerNode, initValue);
   else sender01 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (sender10)  sender10->getData().resize(iNodeSetSender10.size()*iCellSize*sendDataPerNode, initValue);
   else sender10 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (sender11)   sender11->getData().resize(iNodeSetSender11.size()*iCellSize*sendDataPerNode, initValue);
   else sender11 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());

   if (!receiver00) receiver00 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiver01)  receiver01 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiver10)  receiver10 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   if (!receiver11)   receiver11 = VectorTransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::findCFCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes)
{
   int ix1, ix2, ix3;
   LBMReal x1off, x2off, x3off;

   DistributionArray3DPtr  fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();
   BCArray3DPtr bcArray = block.lock()->getKernel()->getBCProcessor()->getBCArray();

   for (ix3 = lMinX3; ix3<=lMaxX3; ix3++)
   {
      for (ix2 = lMinX2; ix2<=lMaxX2; ix2++)
      {
         for (ix1 = lMinX1; ix1<=lMaxX1; ix1++)
         {
            D3Q27ICell icellC;

            int howManySolids = iprocessor->iCellHowManySolids(bcArray, ix1, ix2, ix3);

            if (howManySolids == 0 || howManySolids == 8)
            {
               x1off = 0.0;
               x2off = 0.0;
               x3off = 0.0;
            }
            else
            {
               if (!iprocessor->findNeighborICell(bcArray, fFrom, icellC, bMaxX1, bMaxX2, bMaxX3, ix1, ix2, ix3, x1off, x2off, x3off))
               {
                  std::string err = "For "+block.lock()->toString()+" x1="+UbSystem::toString(ix1)+", x2=" + UbSystem::toString(ix2)+", x3=" + UbSystem::toString(ix3)+
                     " interpolation is not implemented for other direction"+
                     " by using in: "+(std::string)typeid(*this).name()+
                     " or maybe you have a solid on the block boundary";
                  UB_THROW(UbException(UB_EXARGS, err));
               }
            }

            INodeVector inv;
            inv.push_back(ix1 + (int)x1off);
            inv.push_back(ix2 + (int)x2off);
            inv.push_back(ix3 + (int)x3off);
            inv.push_back((int)x1off);
            inv.push_back((int)x2off);
            inv.push_back((int)x3off);
            //inodes.insert(inv);
            inodes.push_back(inv);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
//template< typename VectorTransmitter >
void CoarseToFineNodeSetBlock3DConnector::findCFCells()
{
   using namespace D3Q27System;

   int lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3;

   switch (sendDir)
   {
      //faces
   case E: case W:
      if (sendDir == E)
      {
         lMinX1 = maxX1 - 2;
         lMaxX1 = lMinX1;
      }
      else if (sendDir == W)
      {
         lMinX1 = 1;
         lMaxX1 = lMinX1;
      }

      if (sender00)
      {
         lMinX2 = minX2;
         lMaxX2 = maxHalfX2;
         lMinX3 = minX3;
         lMaxX3 = maxHalfX3;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX2 = minHalfX2;
         lMaxX2 = maxX2 - 1;
         lMinX3 = minX3;
         lMaxX3 = maxHalfX3;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      if (sender01)
      {
         lMinX2 = minX2;
         lMaxX2 = maxHalfX2;
         lMinX3 = minHalfX3;
         lMaxX3 = maxX3 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender01);
      }
      if (sender11)
      {
         lMinX2 = minHalfX2;
         lMaxX2 = maxX2 - 1;
         lMinX3 = minHalfX3;
         lMaxX3 = maxX3 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender11);
      }
      break;
   case N: case S:
      if (sendDir == N)
      {
         lMinX2 = maxX2 - 2;
         lMaxX2 = lMinX2;
      }
      else if (sendDir == S)
      {
         lMinX2 = 1;
         lMaxX2 = lMinX2;
      }

      if (sender00)
      {
         lMinX1 = minX1;
         lMaxX1 = maxHalfX1;
         lMinX3 = minX3;
         lMaxX3 = maxHalfX3;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX1 = minHalfX1;
         lMaxX1 = maxX1 - 1;
         lMinX3 = minX3;
         lMaxX3 = maxHalfX3;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      if (sender01)
      {
         lMinX1 = minX1;
         lMaxX1 = maxHalfX1;
         lMinX3 = minHalfX3;
         lMaxX3 = maxX3 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender01);
      }
      if (sender11)
      {
         lMinX1 = minHalfX1;
         lMaxX1 = maxX1 - 1;
         lMinX3 = minHalfX3;
         lMaxX3 = maxX3 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender11);
      }
      break;
   case T: case B:
      if (sendDir == T)
      {
         lMinX3 = maxX3 - 2;
         lMaxX3 = lMinX3;
      }
      else if (sendDir == B)
      {
         lMinX3 = 1;
         lMaxX3 = lMinX3;
      }

      if (sender00)
      {
         lMinX1 = minX1;
         lMaxX1 = maxHalfX1;
         lMinX2 = minX2;
         lMaxX2 = maxHalfX2;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX1 = minHalfX1;
         lMaxX1 = maxX1 - 1;
         lMinX2 = minX2;
         lMaxX2 = maxHalfX2;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      if (sender01)
      {
         lMinX1 = minX1;
         lMaxX1 = maxHalfX1;
         lMinX2 = minHalfX2;
         lMaxX2 = maxX2 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender01);
      }
      if (sender11)
      {
         lMinX1 = minHalfX1;
         lMaxX1 = maxX1 - 1;
         lMinX2 = minHalfX2;
         lMaxX2 = maxX2 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender11);
      }
      break;
      //edges
      //N-S-E-W
   case NE: case SW: case SE: case NW:
      if (sendDir == NE)
      {
         lMinX1 = maxX1 - 2;
         lMaxX1 = lMinX1 + 1;
         lMinX2 = maxX2 - 2;
         lMaxX2 = lMinX2 + 1;
      }
      else if (sendDir == SW)
      {
         lMinX1 = 0;
         lMaxX1 = lMinX1 + 1;
         lMinX2 = 0;
         lMaxX2 = lMinX2 + 1;
      }
      else if (sendDir == SE)
      {
         lMinX1 = maxX1 - 2;
         lMaxX1 = lMinX1 + 1;
         lMinX2 = 0;
         lMaxX2 = lMinX2 + 1;
      }
      else if (sendDir == NW)
      {
         lMinX1 = 0;
         lMaxX1 = lMinX1 + 1;
         lMinX2 = maxX2 - 2;
         lMaxX2 = lMinX2 + 1;
      }

      if (sender00)
      {
         lMinX3 = minX3;
         lMaxX3 = maxHalfX3;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX3 = minHalfX3;
         lMaxX3 = maxX3 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      break;
      //T-B-E-W
   case TE: case BW: case BE: case TW:
      if (sendDir == TE)
      {
         lMinX1 = maxX1 - 2;
         lMaxX1 = lMinX1 + 1;
         lMinX3 = maxX3 - 2;
         lMaxX3 = lMinX3 + 1;
      }
      else if (sendDir == BW)
      {
         lMinX1 = 0;
         lMaxX1 = lMinX1 + 2;
         lMinX3 = 0;
         lMaxX3 = lMinX3 + 2;
      }
      else if (sendDir == BE)
      {
         lMinX1 = maxX1 - 2;
         lMaxX1 = lMinX1 + 1;
         lMinX3 = 0;
         lMaxX3 = lMinX3 + 1;
      }
      else if (sendDir == TW)
      {
         lMinX1 = 0;
         lMaxX1 = lMinX1 + 1;
         lMinX3 = maxX3 - 2;
         lMaxX3 = lMinX3 + 1;
      }

      if (sender00)
      {
         lMinX2 = minX2;
         lMaxX2 = maxHalfX2;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX2 = minHalfX2;
         lMaxX2 = maxX2 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      break;
      //T-B-N-S
   case TN: case BS: case BN: case TS:
      if (sendDir == TN)
      {
         lMinX2 = maxX2 - 2;
         lMaxX2 = lMinX2 + 1;
         lMinX3 = maxX3 - 3;
         lMaxX3 = lMinX3 + 1;
      }
      else if (sendDir == BS)
      {
         lMinX2 = 0;
         lMaxX2 = lMinX2 + 1;
         lMinX3 = 0;
         lMaxX3 = lMinX3 + 1;
      }
      else if (sendDir == BN)
      {
         lMinX2 = maxX2 - 2;
         lMaxX2 = lMinX2 + 1;
         lMinX3 = 0;
         lMaxX3 = lMinX3 + 1;
      }
      else if (sendDir == TS)
      {
         lMinX2 = 0;
         lMaxX2 = lMinX2 + 1;
         lMinX3 = maxX3 - 2;
         lMaxX3 = lMinX3 + 1;
      }

      if (sender00)
      {
         lMinX1 = minX1;
         lMaxX1 = maxHalfX1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      if (sender10)
      {
         lMinX1 = minHalfX1;
         lMaxX1 = maxX1 - 1;
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender10);
      }
      break;
      //corners
   case TNE: case TNW: case TSE: case TSW: case BNE: case BNW: case BSE: case BSW:
      if (sendDir == TNE)
      {
         lMinX1 = maxX1-2;
         lMaxX1 = maxX1-1;
         lMinX2 = maxX2-2;
         lMaxX2 = maxX2-1;
         lMinX3 = maxX3-2;
         lMaxX3 = maxX3-1;
      }
      else if (sendDir == TNW)
      {
         lMinX1 = 0;
         lMaxX1 = 1;
         lMinX2 = maxX2-2;
         lMaxX2 = maxX2-1;
         lMinX3 = maxX3-2;
         lMaxX3 = maxX3-1;
      }
      else if (sendDir == TSE)
      {
         lMinX1 = maxX1-2;
         lMaxX1 = maxX1-1;
         lMinX2 = 0;
         lMaxX2 = 1;
         lMinX3 = maxX3-2;
         lMaxX3 = maxX3-1;
      }
      else if (sendDir == TSW)
      {
         lMinX1 = 0;
         lMaxX1 = 1;
         lMinX2 = 0;
         lMaxX2 = 1;
         lMinX3 = maxX3-2;
         lMaxX3 = maxX3-1;
      }
      else if (sendDir == BNE)
      {
         lMinX1 = maxX1-2;
         lMaxX1 = maxX1-1;
         lMinX2 = maxX2-2;
         lMaxX2 = maxX2-1;
         lMinX3 = 0;
         lMaxX3 = 1;
      }
      else if (sendDir == BNW)
      {
         lMinX1 = 0;
         lMaxX1 = 1;
         lMinX2 = maxX2-2;
         lMaxX2 = maxX2-1;
         lMinX3 = 0;
         lMaxX3 = 1;
      }
      else if (sendDir == BSE)
      {
         lMinX1 = maxX1-2;
         lMaxX1 = maxX1-1;
         lMinX2 = 0;
         lMaxX2 = 1;
         lMinX3 = 0;
         lMaxX3 = 1;
      }
      else if (sendDir == BSW)
      {
         lMinX1 = 0;
         lMaxX1 = 1;
         lMinX2 = 0;
         lMaxX2 = 1;
         lMinX3 = 0;
         lMaxX3 = 1;
      }
      if (sender00)
      {
         findCFCells(lMinX1, lMinX2, lMinX3, lMaxX1, lMaxX2, lMaxX3, iNodeSetSender00);
      }
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
//template< typename VectorTransmitter >
void CoarseToFineNodeSetBlock3DConnector::fillSendVectors()
{
   using namespace D3Q27System;

   DistributionArray3DPtr  fFrom = block.lock()->getKernel()->getDataSet()->getFdistributions();

   int index00 = 0;
   int index01 = 0;
   int index10 = 0;
   int index11 = 0;

   vector_type& data00 = this->sender00->getData();
   vector_type& data01 = this->sender01->getData();
   vector_type& data10 = this->sender10->getData();
   vector_type& data11 = this->sender11->getData();

   for(INodeVector inode : iNodeSetSender00)
   {
      D3Q27ICell icellC;
      D3Q27ICell icellF;
      iprocessor->readICell(fFrom, icellC, inode[0], inode[1], inode[2]);
      iprocessor->interpolateCoarseToFine(icellC, icellF, inode[3], inode[4], inode[5]);
      writeICellFtoData(data00, index00, icellF);
   }
   for(INodeVector inode : iNodeSetSender01)
   {
      D3Q27ICell icellC;
      D3Q27ICell icellF;
      iprocessor->readICell(fFrom, icellC, inode[0], inode[1], inode[2]);
      iprocessor->interpolateCoarseToFine(icellC, icellF, inode[3], inode[4], inode[5]);
      writeICellFtoData(data01, index01, icellF);
   }
   for(INodeVector inode : iNodeSetSender10)
   {
      D3Q27ICell icellC;
      D3Q27ICell icellF;
      iprocessor->readICell(fFrom, icellC, inode[0], inode[1], inode[2]);
      iprocessor->interpolateCoarseToFine(icellC, icellF, inode[3], inode[4], inode[5]);
      writeICellFtoData(data10, index10, icellF);
   }
   for(INodeVector inode : iNodeSetSender11)
   {
      D3Q27ICell icellC;
      D3Q27ICell icellF;
      iprocessor->readICell(fFrom, icellC, inode[0], inode[1], inode[2]);
      iprocessor->interpolateCoarseToFine(icellC, icellF, inode[3], inode[4], inode[5]);
      writeICellFtoData(data11, index11, icellF);
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::writeICellFtoData(vector_type& data, int& index, D3Q27ICell& icellF)
{
   writeNodeToVector(data, index, icellF.BSW);
   writeNodeToVector(data, index, icellF.BSE);
   writeNodeToVector(data, index, icellF.BNW);
   writeNodeToVector(data, index, icellF.BNE);
   writeNodeToVector(data, index, icellF.TSW);
   writeNodeToVector(data, index, icellF.TSE);
   writeNodeToVector(data, index, icellF.TNW);
   writeNodeToVector(data, index, icellF.TNE);
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::writeNodeToVector(vector_type& data, int& index, LBMReal* inode)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      data[index++] = inode[i];
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::findFCCells(int lMinX1, int lMinX2, int lMinX3, int lMaxX1, int lMaxX2, int lMaxX3, INodeSet &inodes)
{
   int ix1, ix2, ix3;

   for (ix3 = lMinX3; ix3<=lMaxX3; ix3++)
   {
      for (ix2 = lMinX2; ix2<=lMaxX2; ix2++)
      {
         for (ix1 = lMinX1; ix1<=lMaxX1; ix1++)
         {
            INodeVector inv;
            inv.push_back(ix1);
            inv.push_back(ix2);
            inv.push_back(ix3);
            //inodes.insert(inv);
            inodes.push_back(inv);
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
//template< typename VectorTransmitter >
void CoarseToFineNodeSetBlock3DConnector::findFCCells()
{
   using namespace D3Q27System;

   int lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3;
   int lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3;
   int lMin3X1, lMin3X2, lMin3X3, lMax3X1, lMax3X2, lMax3X3;
   int dummy;

   switch (sendDir)
   {

      //////////////////////////////////////////////////////
      //Debug
      //////////////////////////////////////////////////////
      if (block.lock()->getGlobalID() == 2234)
      {
         int test = 0;
      }

      //faces
   case E: case W:
      if (sendDir == E)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = lMin1X1;
      }
      else if (sendDir == W)
      {
         lMin1X1 = 3;
         lMax1X1 = lMin1X1;
      }

      //int TminX1 = lMinX1; int TminX2 = lMinX2; int TminX3 = lMinX3; int TmaxX1 = lMaxX1; int TmaxX2 = lMaxX2; int TmaxX3 = lMaxX3;

      //if (block.lock()->hasInterpolationFlagCF(E))
      //{
      //   if (maxX1==TmaxX1) maxX1 -= 2;
      //}
      //if (block.lock()->hasInterpolationFlagCF(W))
      //{
      //   if (minX1==TminX1) minX1 += 2;
      //}
      //if (block.lock()->hasInterpolationFlagCF(N))
      //{
      //   if (maxX2==TmaxX2)  maxX2 -= 2;
      //}
      //if (block.lock()->hasInterpolationFlagCF(S))
      //{
      //   if (minX2==TminX2)  minX2 += 2;
      //}
      //if (block.lock()->hasInterpolationFlagCF(T))
      //{
      //   if (maxX3==TmaxX3)  maxX3 -= 2;
      //}
      //if (block.lock()->hasInterpolationFlagCF(B))
      //{
      //   if (minX3==TminX3)  minX3 += 2;
      //}
      if (receiver00)
      {
         lMin1X2 = minX2;
         lMax1X2 = maxHalfX2;
         lMin1X3 = minX3;
         lMax1X3 = maxHalfX3;
         getLocalMinMax(dummy, lMin1X2, lMin1X3, dummy, dummy, dummy);
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X2 = minHalfX2;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = minX3;
         lMax1X3 = maxHalfX3;
         getLocalMinMax(dummy, dummy, lMin1X3, dummy, lMax1X2, dummy);
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }
      if (receiver01)
      {
         lMin1X2 = minX2;
         lMax1X2 = maxHalfX2;
         lMin1X3 = minHalfX3;
         lMax1X3 = maxX3 - 1;
         getLocalMinMax(dummy, lMin1X2, dummy, dummy, dummy, lMax1X3);
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver01);
      }
      if (receiver11)
      {
         lMin1X2 = minHalfX2;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = minHalfX3;
         lMax1X3 = maxX3 - 1;
         getLocalMinMax(dummy, dummy, dummy, dummy, lMax1X2, lMax1X3);
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver11);
      }
      break;
   case N: case S:
      if (sendDir == N)
      {
         lMin1X2 = maxX2 - 3;
         lMax1X2 = lMin1X2;
      }
      else if (sendDir == S)
      {
         lMin1X2 = 3;
         lMax1X2 = lMin1X2;
      }

      if (receiver00)
      {
         lMin1X1 = minX1;
         lMax1X1 = maxHalfX1;
         lMin1X3 = minX3;
         lMax1X3 = maxHalfX3;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X1 = minHalfX1;
         lMax1X1 = maxX1 - 1;
         lMin1X3 = minX3;
         lMax1X3 = maxHalfX3;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }
      if (receiver01)
      {
         lMin1X1 = minX1;
         lMax1X1 = maxHalfX1;
         lMin1X3 = minHalfX3;
         lMax1X3 = maxX3 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver01);
      }
      if (receiver11)
      {
         lMin1X1 = minHalfX1;
         lMax1X1 = maxX1 - 1;
         lMin1X3 = minHalfX3;
         lMax1X3 = maxX3 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver11);
      }
      break;
   case T: case B:
      if (sendDir == T)
      {
         lMin1X3 = maxX3 - 3;
         lMax1X3 = lMin1X3;
      }
      else if (sendDir == B)
      {
         lMin1X3 = 3;
         lMax1X3 = lMin1X3;
      }

      if (receiver00)
      {
         lMin1X1 = minX1;
         lMax1X1 = maxHalfX1;
         lMin1X2 = minX2;
         lMax1X2 = maxHalfX2;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X1 = minHalfX1;
         lMax1X1 = maxX1 - 1;
         lMin1X2 = minX2;
         lMax1X2 = maxHalfX2;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }
      if (receiver01)
      {
         lMin1X1 = minX1;
         lMax1X1 = maxHalfX1;
         lMin1X2 = minHalfX2;
         lMax1X2 = maxX2 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver01);
      }
      if (receiver11)
      {
         lMin1X1 = minHalfX1;
         lMax1X1 = maxX1 - 1;
         lMin1X2 = minHalfX2;
         lMax1X2 = maxX2 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver11);
      }
      break;
      //edges
      //N-S-E-W
   case NE: case SW: case SE: case NW:
      if (sendDir == NE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = lMin1X1 + 2;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = lMin1X2;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = lMin2X1;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = lMin2X2 + 2;
      }
      else if (sendDir == SW)
      {
         lMin1X1 = 1;
         lMax1X1 = lMin1X1 + 2;
         lMin1X2 = 3;
         lMax1X2 = lMin1X2;

         lMin2X1 = 3;
         lMax2X1 = lMin2X1;
         lMin2X2 = 1;
         lMax2X2 = lMin2X2 + 2;
      }
      else if (sendDir == SE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = lMin1X1 + 2;
         lMin1X2 = 3;
         lMax1X2 = lMin1X2;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = lMin2X1;
         lMin2X2 = 1;
         lMax2X2 = lMin2X2 + 2;
      }
      else if (sendDir == NW)
      {
         lMin1X1 = 1;
         lMax1X1 = lMin1X1 + 2;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = lMin1X2;

         lMin2X1 = 3;
         lMax2X1 = lMin2X1;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = lMin2X2 + 2;
      }

      if (receiver00)
      {
         lMin1X3 = minX3;
         lMax1X3 = maxHalfX3;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X3 = minHalfX3;
         lMax1X3 = maxX3 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }

      if (receiver00)
      {
         lMin2X3 = minX3;
         lMax2X3 = maxHalfX3;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin2X3 = minHalfX3;
         lMax2X3 = maxX3 - 1;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver10);
      }
      break;
      //T-B-E-W
   case TE: case BW: case BE: case TW:
      if (sendDir == TE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = lMin1X1 + 2;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = lMin1X3;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = lMin2X1;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == BW)
      {
         lMin1X1 = 1;
         lMax1X1 = lMin1X1 + 2;
         lMin1X3 = 3;
         lMax1X3 = lMin1X3;

         lMin2X1 = 3;
         lMax2X1 = lMin2X1;
         lMin2X3 = 1;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == BE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = lMin1X1 + 2;
         lMin1X3 = 3;
         lMax1X3 = lMin1X3;


         lMin2X1 = maxX1 - 3;
         lMax2X1 = lMin2X1;
         lMin2X3 = 1;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == TW)
      {
         lMin1X1 = 1;
         lMax1X1 = lMin1X1 + 2;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = lMin1X3;

         lMin2X1 = 3;
         lMax2X1 = lMin2X1;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = lMin2X3 + 2;
      }

      if (receiver00)
      {
         lMin1X2 = minX2;
         lMax1X2 = maxHalfX2;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X2 = minHalfX2;
         lMax1X2 = maxX2 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }

      if (receiver00)
      {
         lMin2X2 = minX2;
         lMax2X2 = maxHalfX2;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin2X2 = minHalfX2;
         lMax2X2 = maxX2 - 1;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver10);
      }
      break;
      //T-B-N-S
   case TN: case BS: case BN: case TS:
      if (sendDir == TN)
      {
         lMin1X2 = maxX2 - 3;
         lMax1X2 = lMin1X2 + 2;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = lMin1X3;

         lMin2X2 = maxX2 - 3;
         lMax2X2 = lMin2X2;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == BS)
      {
         lMin1X2 = 1;
         lMax1X2 = lMin1X2 + 2;
         lMin1X3 = 3;
         lMax1X3 = lMin1X3;

         lMin2X2 = 3;
         lMax2X2 = lMin2X2;
         lMin2X3 = 1;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == BN)
      {
         lMin1X2 = maxX2 - 3;
         lMax1X2 = lMin1X2 + 2;
         lMin1X3 = 3;
         lMax1X3 = lMin1X3;

         lMin2X2 = maxX2 - 3;
         lMax2X2 = lMin2X2;
         lMin2X3 = 1;
         lMax2X3 = lMin2X3 + 2;
      }
      else if (sendDir == TS)
      {
         lMin1X2 = 1;
         lMax1X2 = lMin1X2 + 2;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = lMin1X3;

         lMin2X2 = 3;
         lMax2X2 = lMin2X2;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = lMin2X3 + 2;
      }

      if (receiver00)
      {
         lMin1X1 = minX1;
         lMax1X1 = maxHalfX1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin1X1 = minHalfX1;
         lMax1X1 = maxX1 - 1;
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver10);
      }

      if (receiver00)
      {
         lMin2X1 = minX1;
         lMax2X1 = maxHalfX1;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver00);
      }
      if (receiver10)
      {
         lMin2X1 = minHalfX1;
         lMax2X1 = maxX1 - 1;
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver10);
      }
      break;
      //corners
   case TNE: case TNW: case TSE: case TSW: case BNE: case BNW: case BSE: case BSW:
      if (sendDir == TNE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = maxX1 - 2;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = maxX3 - 1;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = maxX1 - 1;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = maxX2 - 2;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = maxX3 - 1;

         lMin3X1 = maxX1 - 3;
         lMax3X1 = maxX1 - 1;
         lMin3X2 = maxX2 - 3;
         lMax3X2 = maxX2 - 1;
         lMin3X3 = maxX3 - 3;
         lMax3X3 = maxX3 - 2;
      }
      else if (sendDir == TNW)
      {
         lMin1X1 = 3;
         lMax1X1 = 3;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = maxX3 - 1;

         lMin2X1 = 1;
         lMax2X1 = 3;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = maxX2 - 2;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = maxX3;

         lMin3X1 = 1;
         lMax3X1 = 3;
         lMin3X2 = maxX2 - 3;
         lMax3X2 = maxX2 - 1;
         lMin3X3 = maxX3 - 3;
         lMax3X3 = maxX3 - 2;
      }
      else if (sendDir == TSE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = maxX1 - 2;
         lMin1X2 = 1;
         lMax1X2 = 3;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = maxX3;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = maxX1 - 1;
         lMin2X2 = 3;
         lMax2X2 = 3;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = maxX3;

         lMin3X1 = maxX1 - 3;
         lMax3X1 = maxX1 - 1;
         lMin3X2 = 1;
         lMax3X2 = 3;
         lMin3X3 = maxX3 - 3;
         lMax3X3 = maxX3 - 2;
      }
      else if (sendDir == TSW)
      {
         lMin1X1 = 3;
         lMax1X1 = 3;
         lMin1X2 = 1;
         lMax1X2 = 3;
         lMin1X3 = maxX3 - 3;
         lMax1X3 = maxX3 - 1;

         lMin2X1 = 1;
         lMax2X1 = 3;
         lMin2X2 = 3;
         lMax2X2 = 3;
         lMin2X3 = maxX3 - 3;
         lMax2X3 = maxX3 - 1;

         lMin3X1 = 1;
         lMax3X1 = 3;
         lMin3X2 = 1;
         lMax3X2 = 3;
         lMin3X3 = maxX3 - 3;
         lMax3X3 = maxX3 - 2;
      }
      else if (sendDir == BNE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = maxX1 - 2;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = 1;
         lMax1X3 = 3;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = maxX1 - 1;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = maxX2 - 2;
         lMin2X3 = 1;
         lMax2X3 = 3;

         lMin3X1 = maxX1 - 3;
         lMax3X1 = maxX1 - 1;
         lMin3X2 = maxX2 - 3;
         lMax3X2 = maxX2 - 1;
         lMin3X3 = 3;
         lMax3X3 = 3;
      }
      else if (sendDir == BNW)
      {
         lMin1X1 = 3;
         lMax1X1 = 3;
         lMin1X2 = maxX2 - 3;
         lMax1X2 = maxX2 - 1;
         lMin1X3 = 1;
         lMax1X3 = 3;

         lMin2X1 = 1;
         lMax2X1 = 3;
         lMin2X2 = maxX2 - 3;
         lMax2X2 = maxX2 - 2;
         lMin2X3 = 1;
         lMax2X3 = 3;

         lMin3X1 = 1;
         lMax3X1 = 3;
         lMin3X2 = maxX2 - 3;
         lMax3X2 = maxX2 - 1;
         lMin3X3 = 3;
         lMax3X3 = 3;
      }
      else if (sendDir == BSE)
      {
         lMin1X1 = maxX1 - 3;
         lMax1X1 = maxX1 - 2;
         lMin1X2 = 1;
         lMax1X2 = 3;
         lMin1X3 = 1;
         lMax1X3 = 3;

         lMin2X1 = maxX1 - 3;
         lMax2X1 = maxX1 - 1;
         lMin2X2 = 3;
         lMax2X2 = 3;
         lMin2X3 = 1;
         lMax2X3 = 3;

         lMin3X1 = maxX1 - 3;
         lMax3X1 = maxX1 - 1;
         lMin3X2 = 1;
         lMax3X2 = 3;
         lMin3X3 = 3;
         lMax3X3 = 3;
      }
      else if (sendDir == BSW)
      {
         lMin1X1 = 3;
         lMax1X1 = 3;
         lMin1X2 = 1;
         lMax1X2 = 3;
         lMin1X3 = 1;
         lMax1X3 = 3;

         lMin2X1 = 1;
         lMax2X1 = 3;
         lMin2X2 = 3;
         lMax2X2 = 3;
         lMin2X3 = 1;
         lMax2X3 = 3;

         lMin3X1 = 1;
         lMax3X1 = 3;
         lMin3X2 = 1;
         lMax3X2 = 3;
         lMin3X3 = 3;
         lMax3X3 = 3;
      }
      if (receiver00)
      {
         findFCCells(lMin1X1, lMin1X2, lMin1X3, lMax1X1, lMax1X2, lMax1X3, iNodeSetReceiver00);
      }
      if (receiver00)
      {
         findFCCells(lMin2X1, lMin2X2, lMin2X3, lMax2X1, lMax2X2, lMax2X3, iNodeSetReceiver00);
      }
      if (receiver00)
      {
         findFCCells(lMin3X1, lMin3X2, lMin3X3, lMax3X1, lMax3X2, lMax3X3, iNodeSetReceiver00);
      }
      break;
   }
}
//////////////////////////////////////////////////////////////////////////
//template< typename VectorTransmitter >
void CoarseToFineNodeSetBlock3DConnector::distributeReceiveVectors()
{
   using namespace D3Q27System;

   DistributionArray3DPtr  fTo = block.lock()->getKernel()->getDataSet()->getFdistributions();

   int index00 = 0;
   int index01 = 0;
   int index10 = 0;
   int index11 = 0;

   vector_type& data00 = this->receiver00->getData();
   vector_type& data01 = this->receiver01->getData();
   vector_type& data10 = this->receiver10->getData();
   vector_type& data11 = this->receiver11->getData();

   for(INodeVector inode : iNodeSetReceiver00)
   {
      LBMReal icellC[27];
      this->readICellCfromData(data00, index00, icellC);
      iprocessor->writeINodeInv(fTo, icellC, inode[0], inode[1], inode[2]);
   }
   for(INodeVector inode : iNodeSetReceiver01)
   {
      LBMReal icellC[27];
      this->readICellCfromData(data01, index01, icellC);
      iprocessor->writeINodeInv(fTo, icellC, inode[0], inode[1], inode[2]);
   }
   for(INodeVector inode : iNodeSetReceiver10)
   {
      LBMReal icellC[27];
      this->readICellCfromData(data10, index10, icellC);
      iprocessor->writeINodeInv(fTo, icellC, inode[0], inode[1], inode[2]);
   }
   for(INodeVector inode : iNodeSetReceiver11)
   {
      LBMReal icellC[27];
      this->readICellCfromData(data11, index11, icellC);
      iprocessor->writeINodeInv(fTo, icellC, inode[0], inode[1], inode[2]);
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::readICellCfromData(vector_type& data, int& index, LBMReal* icellC)
{
   for (int i = D3Q27System::STARTF; i < D3Q27System::ENDF+1; i++)
   {
      icellC[i] = data[index++];
   }
}
//////////////////////////////////////////////////////////////////////////

void CoarseToFineNodeSetBlock3DConnector::getLocalMinMax(int& minX1, int& minX2, int& minX3, int& maxX1, int& maxX2, int& maxX3)
{
   using namespace D3Q27System;
   int TminX1 = minX1; int TminX2 = minX2; int TminX3 = minX3; int TmaxX1 = maxX1; int TmaxX2 = maxX2; int TmaxX3 = maxX3;

   if (block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1==TmaxX1) maxX1 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1==TminX1) minX1 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(N))
   {
      if (maxX2==TmaxX2)  maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX2==TminX2)  minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX3==TmaxX3)  maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX3==TminX3)  minX3 += 2;
   }

   //E-W-N-S
   if (block.lock()->hasInterpolationFlagCF(NE) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(E))
   {
      if (maxX1==TmaxX1) maxX1 -= 2;
      if (maxX2==TmaxX2) maxX2 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (minX1==TminX1) minX1 += 2;
      if (minX2==TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(SE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(S))
   {
      if (maxX1==TmaxX1) maxX1 -= 2;
      if (minX2==TminX2) minX2 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(NW) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(W))
   {
      if (minX1==TminX1) minX1 += 2;
      if (maxX2==TmaxX2) maxX2 -= 2;
   }

   //	////T-B-E-W
   if (block.lock()->hasInterpolationFlagCF(TE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX1==TmaxX1) maxX1 -= 2;
      if (maxX3==TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX1==TminX1) minX1 += 2;
      if (minX3==TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BE) && !block.lock()->hasInterpolationFlagCF(E) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX1==TmaxX1) maxX1 -= 2;
      if (minX3==TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TW) && !block.lock()->hasInterpolationFlagCF(W) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX1==TminX1) minX1 += 2;
      if (maxX3==TmaxX3) maxX3 -= 2;
   }


   ////T-B-N-S
   if (block.lock()->hasInterpolationFlagCF(TN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (maxX2==TmaxX2) maxX2 -= 2;
      if (maxX3==TmaxX3) maxX3 -= 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (minX2==TminX2) minX2 += 2;
      if (minX3==TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(BN) && !block.lock()->hasInterpolationFlagCF(N) && !block.lock()->hasInterpolationFlagCF(B))
   {
      if (maxX2==TmaxX2) maxX2 -= 2;
      if (minX3==TminX3) minX3 += 2;
   }
   if (block.lock()->hasInterpolationFlagCF(TS) && !block.lock()->hasInterpolationFlagCF(S) && !block.lock()->hasInterpolationFlagCF(T))
   {
      if (minX2==TminX2) minX2 += 2;
      if (maxX3==TmaxX3) maxX3 -= 2;
   }
}


