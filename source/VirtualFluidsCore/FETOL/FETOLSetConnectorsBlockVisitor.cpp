#if defined VF_MPI && defined VF_FETOL

#include "FETOLSetConnectorsBlockVisitor.h"
#include "D3Q27ETFullVectorConnector.h"
#include "D3Q27ETFullDirectConnector.h"
#include "D3Q27ETFullDirectConnector2.h"
#include "D3Q27ETFullVectorConnector.h"
#include "D3Q27ETCFVectorConnector.h"
#include "D3Q27ETFCVectorConnector.h"
#include "D3Q27ETCFOffVectorConnector.h"
#include "D3Q27ETFCOffVectorConnector.h"
#include "Grid3DSystem.h"
#include "D3Q27ETDirectConnector.h"
#include <basics/transmitter/TbTransmitterLocal.h>

FETOLSetConnectorsBlockVisitor::FETOLSetConnectorsBlockVisitor(CommunicatorPtr comm, bool fullConnector, int dirs, 
                                                               LBMReal nue, D3Q27InterpolationProcessorPtr iProcessor) :
                                                               Block3DVisitor(0, Grid3DSystem::MAXLEVEL), 
                                                               fullConnector(fullConnector),
                                                               dirs(dirs),
                                                               nue(nue),
                                                               iProcessor(iProcessor)
{
   this->comm = boost::dynamic_pointer_cast<FETOLCommunicator>(comm);
   gridBundle = this->comm->getBundleID();
   gridRank = comm->getProcessID();
   UBLOG(logDEBUG5, "D3Q27BondSetConnectorsBlockVisitor: gridBundle = "<<gridBundle<<" gridRank = "<<gridRank);
}
//////////////////////////////////////////////////////////////////////////
FETOLSetConnectorsBlockVisitor::~FETOLSetConnectorsBlockVisitor(void)
{
}
//////////////////////////////////////////////////////////////////////////
void FETOLSetConnectorsBlockVisitor::visit(Grid3DPtr grid, Block3DPtr block)
{
   if(!block) return;

   //UBLOG(logDEBUG5, "D3Q27BondSetConnectorsBlockVisitor::visit() - start");
   //gridBundle = comm->getBundleID();
   //grid->setBundle(gridBundle);

   //gridRank = comm->getProcessID();
   //grid->setRank(gridRank);
   
   setSameLevelConnectors(grid, block);

   //if(grid->getFinestInitializedLevel() > grid->getCoarsestInitializedLevel())
   //   setInterpolationConnectors(grid, block);

   //UBLOG(logDEBUG5, "D3Q27BondSetConnectorsBlockVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void FETOLSetConnectorsBlockVisitor::setSameLevelConnectors(Grid3DPtr grid, Block3DPtr block)
{
   int blockBundle = block->getBundle();
   int blockRank = block->getRank();
   if (gridBundle == blockBundle && gridRank == blockRank && block->isActive())
   {
      block->clearWeight();
      std::vector<Block3DPtr> neighbors; 
      int ix1 = block->getX1();
      int ix2 = block->getX2();
      int ix3 = block->getX3();
      int level = block->getLevel();

      for( int dir = 0; dir < dirs; dir++)
      {
         Block3DPtr neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);

         if(neighBlock)
         {
            int neighBlockBundle = neighBlock->getBundle();
            int neighBlockRank = neighBlock->getRank();
            if(blockBundle == neighBlockBundle && blockRank == neighBlockRank && neighBlock->isActive())
            {
               Block3DConnectorPtr connector;
               connector = Block3DConnectorPtr(new D3Q27ETFullDirectConnector2( block, neighBlock, dir));
               block->setConnector(connector);
            }
            else if(blockBundle == neighBlockBundle && blockRank != neighBlockRank && neighBlock->isActive())
            {
               setRemoteConnectors(block, neighBlock, dir); 
            }
            else if(blockBundle != neighBlockBundle && neighBlock->isActive())
            {
               setBundleConnectors(block, neighBlock, dir);  
            }
         }
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void FETOLSetConnectorsBlockVisitor::setRemoteConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir)
{
   CreateTransmittersHelper helper;
   CreateTransmittersHelper::TransmitterPtr sender, receiver;
   helper.createTransmitters(sblock, tblock, dir, CreateTransmittersHelper::NONE, sender, receiver, comm, CreateTransmittersHelper::MPI2BOND);

   Block3DConnectorPtr connector;
   connector = Block3DConnectorPtr(new D3Q27ETFullVectorConnector< TbTransmitter< CbVector< LBMReal > > >(sblock, sender, receiver, dir));
   connector->setTransmitterType(Block3DConnector::MPI);
   UBLOG(logDEBUG5,"setTransmitterType: "<< Block3DConnector::MPI);
   sblock->setConnector(connector);
   UBLOG(logDEBUG5, "D3Q27BondSetConnectorsBlockVisitor::setRemoteConnectors() - sblock = "<<sblock->toString());
}
//////////////////////////////////////////////////////////////////////////
void FETOLSetConnectorsBlockVisitor::setBundleConnectors(Block3DPtr sblock, Block3DPtr tblock, int dir)
{
   CreateTransmittersHelper helper;
   CreateTransmittersHelper::TransmitterPtr sender, receiver;
   bool MPIpool = true;
   bool orthogonal = false;
   helper.createTransmitters(sblock, tblock, dir, CreateTransmittersHelper::NONE, sender, receiver, comm, CreateTransmittersHelper::BOND);

   Block3DConnectorPtr connector;
   connector = Block3DConnectorPtr(new D3Q27ETFullVectorConnector< TbTransmitter< CbVector< LBMReal > > >(sblock, sender, receiver, dir));
   connector->setTransmitterType(Block3DConnector::BOND);
   UBLOG(logDEBUG5,"setTransmitterType: "<< Block3DConnector::BOND);
   sblock->setConnector(connector);
   UBLOG(logDEBUG5, "D3Q27BondSetConnectorsBlockVisitor::setBundleConnectors() - sblock = "<<sblock->toString());
}
//////////////////////////////////////////////////////////////////////////
//void D3Q27BondSetConnectorsBlockVisitor::setInterpolationConnectors(Grid3DPtr grid, Block3DPtr block)
//{
//   //int blockRank = block->getRank();
//
//   //search for all blocks with different ranks
//   if (block->hasInterpolationFlagCF() && block->isActive())
//   {
//      int fbx1 = block->getX1() << 1;
//      int fbx2 = block->getX2() << 1;
//      int fbx3 = block->getX3() << 1;
//      int level = block->getLevel() + 1;
//
//      if( block->hasInterpolationFlagCF(D3Q27System::E))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1+1,fbx2,fbx3,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::E);
//      }
//      if( block->hasInterpolationFlagCF(D3Q27System::W))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1,fbx2+1,fbx3,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::W);
//      }
//      if( block->hasInterpolationFlagCF(D3Q27System::N))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::N);
//      }
//      if( block->hasInterpolationFlagCF(D3Q27System::S))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::S);
//      }
//      if( block->hasInterpolationFlagCF(D3Q27System::T))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1,fbx2,fbx3+1,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3+1,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3+1,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3+1,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::T);
//      }
//      if( block->hasInterpolationFlagCF(D3Q27System::B))
//      {
//         Block3DPtr fblockSW = grid->getBlock(fbx1,fbx2,fbx3,level);
//         Block3DPtr fblockSE = grid->getBlock(fbx1+1,fbx2,fbx3,level);
//         Block3DPtr fblockNW = grid->getBlock(fbx1,fbx2+1,fbx3,level);
//         Block3DPtr fblockNE = grid->getBlock(fbx1+1,fbx2+1,fbx3,level);
//
//         setInterpolationConnectors(grid, fblockSW, fblockSE, fblockNW, fblockNE, block, D3Q27System::B);
//      }
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void D3Q27BondSetConnectorsBlockVisitor::setInterpolationConnectors(Grid3DPtr grid, Block3DPtr fBlockSW, Block3DPtr fBlockSE, Block3DPtr fBlockNW, Block3DPtr fBlockNE, Block3DPtr cBlock, int dir)
//{
//   int fBlockSWRank = -999, fBlockSERank = -999, fBlockNWRank = -999, fBlockNERank = -999;
//   if(fBlockSW) fBlockSWRank = fBlockSW->getRank();
//   if(fBlockNW) fBlockNWRank = fBlockNW->getRank();
//   if(fBlockSE) fBlockSERank = fBlockSE->getRank();
//   if(fBlockNE) fBlockNERank = fBlockNE->getRank();
//   int cBlockRank   = cBlock->getRank();
//
//   LBMReal omegaF;
//   if(fBlockSW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSW->getLevel());
//   if(fBlockNW) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNW->getLevel());
//   if(fBlockSE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockSE->getLevel());
//   if(fBlockNE) omegaF =LBMSystem::calcCollisionFactor(nue, fBlockNE->getLevel());
//   LBMReal omegaC = LBMSystem::calcCollisionFactor(nue, cBlock->getLevel());
//   iProcessor->setOmegas(omegaC, omegaF);
//
//   D3Q27InterpolationProcessorPtr cIProcessor(iProcessor->clone());
//   D3Q27InterpolationProcessorPtr fIProcessorSW(iProcessor->clone());
//   D3Q27InterpolationProcessorPtr fIProcessorSE(iProcessor->clone());
//   D3Q27InterpolationProcessorPtr fIProcessorNW(iProcessor->clone());
//   D3Q27InterpolationProcessorPtr fIProcessorNE(iProcessor->clone());
//
//   D3Q27CreateTransmittersHelper::TransmitterPtr senderCFevenEvenSW, receiverCFevenEvenSW, 
//      senderCFevenOddNW,  receiverCFevenOddNW, 
//      senderCFoddEvenSE,  receiverCFoddEvenSE, 
//      senderCFoddOddNE,   receiverCFoddOddNE,
//      senderFCevenEvenSW, receiverFCevenEvenSW, 
//      senderFCevenOddNW,  receiverFCevenOddNW, 
//      senderFCoddEvenSE,  receiverFCoddEvenSE, 
//      senderFCoddOddNE,   receiverFCoddOddNE;
//
//   if(fBlockSW) createTransmitters(cBlock, fBlockSW, dir, senderCFevenEvenSW, receiverCFevenEvenSW, senderFCevenEvenSW, receiverFCevenEvenSW);
//   if(fBlockNW) createTransmitters(cBlock, fBlockNW, dir, senderCFevenOddNW, receiverCFevenOddNW, senderFCevenOddNW, receiverFCevenOddNW);
//   if(fBlockSE) createTransmitters(cBlock, fBlockSE, dir, senderCFoddEvenSE, receiverCFoddEvenSE, senderFCoddEvenSE, receiverFCoddEvenSE);
//   if(fBlockNE) createTransmitters(cBlock, fBlockNE, dir, senderCFoddOddNE, receiverCFoddOddNE, senderFCoddOddNE, receiverFCoddOddNE);
//
//   if(cBlockRank == gridRank)
//   {
//      Block3DConnectorPtr connector( new D3Q27ETCFOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(cBlock, 
//         senderCFevenEvenSW, receiverCFevenEvenSW, senderCFevenOddNW,  receiverCFevenOddNW, 
//         senderCFoddEvenSE,  receiverCFoddEvenSE,  senderCFoddOddNE,   receiverCFoddOddNE, 
//         dir, cIProcessor) );
//      cBlock->setConnector(connector);
//   }
//   if(fBlockSW && fBlockSWRank == gridRank)
//   {
//      Block3DConnectorPtr connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockSW, 
//         senderFCevenEvenSW, receiverFCevenEvenSW, dir, fIProcessorSW, EvenEvenSW) );
//      fBlockSW->setConnector(connector);
//   }
//   if(fBlockNW && fBlockNWRank == gridRank)
//   {
//      Block3DConnectorPtr connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockNW, 
//         senderFCevenOddNW, receiverFCevenOddNW, dir, fIProcessorNW, EvenOddNW) );
//      fBlockNW->setConnector(connector);
//   }
//   if(fBlockSE && fBlockSERank == gridRank)
//   {
//      Block3DConnectorPtr connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockSE, 
//         senderFCoddEvenSE, receiverFCoddEvenSE, dir, fIProcessorSE, OddEvenSE) );
//      fBlockSE->setConnector(connector);
//   }
//   if(fBlockNE && fBlockNERank == gridRank)
//   {
//      Block3DConnectorPtr connector( new D3Q27ETFCOffVectorConnector< TbTransmitter< CbVector< LBMReal > > >(fBlockNE, 
//         senderFCoddOddNE, receiverFCoddOddNE, dir, fIProcessorNE, OddOddNE) );
//      fBlockNE->setConnector(connector);
//   }
//}
////////////////////////////////////////////////////////////////////////////
//void D3Q27BondSetConnectorsBlockVisitor::createTransmitters(Block3DPtr cBlock, Block3DPtr fBlock, int dir, 
//                                                         D3Q27CreateTransmittersHelper::TransmitterPtr& senderCF, 
//                                                         D3Q27CreateTransmittersHelper::TransmitterPtr& receiverCF, 
//                                                         D3Q27CreateTransmittersHelper::TransmitterPtr& senderFC, 
//                                                         D3Q27CreateTransmittersHelper::TransmitterPtr& receiverFC)
//{
//   D3Q27CreateTransmittersHelper helper;
//   bool MPIpool = true;
//   bool orthogonal = false;
//   int fBlockRank = fBlock->getRank();
//   int cBlockRank = cBlock->getRank();
//   if(fBlockRank == cBlockRank && fBlockRank == gridRank)
//   {
//      senderCF = receiverFC = D3Q27CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< LBMReal > >());
//      senderFC = receiverCF = D3Q27CreateTransmittersHelper::TransmitterPtr( new TbLocalTransmitter< CbVector< LBMReal > >());
//   }
//   else if(cBlockRank == gridRank)
//   {
//      helper.createTransmitters(cBlock, fBlock, dir, senderCF, receiverCF, comm);
//   }
//   else if(fBlockRank == gridRank)
//   {
//      helper.createTransmitters(fBlock, cBlock, dir, senderFC, receiverFC, comm);
//   }
//}

#endif
