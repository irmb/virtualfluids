#include "ConnectorBlockVisitor.h"
#include "Grid3DSystem.h"
#include "ConnectorFactory.h"
#include "InterpolationProcessor.h"
#include "Communicator.h"
#include "Grid3D.h"

ConnectorBlockVisitor::ConnectorBlockVisitor(SPtr<Communicator> comm, LBMReal nu, InterpolationProcessorPtr iProcessor, SPtr<ConnectorFactory> cFactory) :
   Block3DVisitor(0, Grid3DSystem::MAXLEVEL),
   comm(comm),
   nu(nu),
   iProcessor(iProcessor),
   cFactory(cFactory)
{
}
//////////////////////////////////////////////////////////////////////////
ConnectorBlockVisitor::~ConnectorBlockVisitor(void)
{
}
//////////////////////////////////////////////////////////////////////////
void ConnectorBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (!block) return;

   UBLOG(logDEBUG5, "ConnectorBlockVisitor::visit() - start");
   UBLOG(logDEBUG5, block->toString());

   gridRank = comm->getProcessID();
   grid->setRank(gridRank);

   setSameLevelConnectors(grid, block);

   if (grid->getFinestInitializedLevel() > grid->getCoarsestInitializedLevel())
      setInterpolationConnectors(grid, block);

   if (block->getGlobalID()==2234)
   {
      UBLOG(logINFO, block->toString());
   }

   UBLOG(logDEBUG5, "ConnectorBlockVisitor::visit() - end");
}
//////////////////////////////////////////////////////////////////////////
void ConnectorBlockVisitor::setSameLevelConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (block->getGlobalID()==2234)
   {
      UBLOG(logINFO, block->toString());
   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setSameLevelConnectors() - start");
   int blockRank = block->getRank();
   if (gridRank == blockRank && block->isActive())
   {
      block->clearWeight();
      std::vector<SPtr<Block3D>> neighbors;
      int ix1 = block->getX1();
      int ix2 = block->getX2();
      int ix3 = block->getX3();
      int level = block->getLevel();

      for (int dir = 0; dir < D3Q27System::ENDDIR; dir++)
      {
         SPtr<Block3D> neighBlock = grid->getNeighborBlock(dir, ix1, ix2, ix3, level);

         if (neighBlock)
         {
            int neighBlockRank = neighBlock->getRank();
            if (blockRank == neighBlockRank && neighBlock->isActive())
            {
               SPtr<Block3DConnector> connector;
               connector = cFactory->createSameLevelDirectConnector(block, neighBlock, dir);
               block->setConnector(connector);
            }
            else if (blockRank != neighBlockRank && neighBlock->isActive())
            {
               setRemoteConnectors(block, neighBlock, dir);
            }
         }
      }
   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setSameLevelConnectors() - end");
   if (block->getGlobalID()==2234)
   {
      UBLOG(logINFO, block->toString());
   }
}
//////////////////////////////////////////////////////////////////////////
void ConnectorBlockVisitor::setRemoteConnectors(SPtr<Block3D> sblock, SPtr<Block3D> tblock, int dir)
{
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setRemoteConnectors() - start");
   CreateTransmittersHelper helper;
   CreateTransmittersHelper::TransmitterPtr sender, receiver;
   helper.createTransmitters(sblock, tblock, dir, CreateTransmittersHelper::NONE, sender, receiver, comm, CreateTransmittersHelper::MPI);


   SPtr<Block3DConnector> connector;
   connector = cFactory->createSameLevelVectorConnector(sblock, sender, receiver, dir);
   sblock->setConnector(connector);
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setRemoteConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////
void ConnectorBlockVisitor::setInterpolationConnectors(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
   if (block->getGlobalID()==2234)
   {
      UBLOG(logINFO, block->toString());
   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setInterpolationConnectors() - start");
//   int blockRank = block->getRank();

   //search for all blocks with different ranks
   if (block->hasInterpolationFlagCF() && block->isActive())
   {
      int fbx1 = block->getX1() << 1;
      int fbx2 = block->getX2() << 1;
      int fbx3 = block->getX3() << 1;
      int level = block->getLevel() + 1;

      if (block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::E);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::W);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::N))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::N);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::S);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::T))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::T);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::B))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblock01 = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock11 = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::B);
      }

      //////NE-NW-SE-SW
      if (block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::NE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::SW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::SE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::NW);
      }

      /////////TE-BW-BE-TW 1-0
      if (block->hasInterpolationFlagCF(D3Q27System::TE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::T))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2+0, fbx3+1, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::B))
      {

         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+0, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BE)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::B))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2+0, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TW)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::T))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+0, fbx3+1, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TW);
      }

      //////TN-BS-BN-TS
      if (block->hasInterpolationFlagCF(D3Q27System::TN)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::T))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+0, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TN);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::S) && !block->hasInterpolationFlagCF(D3Q27System::B))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+0, fbx2, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BS);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::B))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+0, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BN);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::S) && !block->hasInterpolationFlagCF(D3Q27System::T))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+0, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock10 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock01;
         SPtr<Block3D> fblock11;

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TS);
      }




      //////corners
      if (block->hasInterpolationFlagCF(D3Q27System::TNE)&&!block->hasInterpolationFlagCF(D3Q27System::TE)&&!block->hasInterpolationFlagCF(D3Q27System::TN)&&!block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::T)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TNE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TSW)&&!block->hasInterpolationFlagCF(D3Q27System::TW)&&!block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TSW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TSE)&&!block->hasInterpolationFlagCF(D3Q27System::TE)&&!block->hasInterpolationFlagCF(D3Q27System::TS)&& !block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TSE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::TNW)&&!block->hasInterpolationFlagCF(D3Q27System::TW)&&!block->hasInterpolationFlagCF(D3Q27System::TN)&& !block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::T)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::TNW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BNE)&&!block->hasInterpolationFlagCF(D3Q27System::BE)&&!block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::NE)&&!block->hasInterpolationFlagCF(D3Q27System::B)&&!block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::E))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1+1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BNE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BSW)&& !block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::SW)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::W) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1, fbx2, fbx3, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BSW);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BSE)&& !block->hasInterpolationFlagCF(D3Q27System::BS)&& !block->hasInterpolationFlagCF(D3Q27System::BE)&& !block->hasInterpolationFlagCF(D3Q27System::SE)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::E) && !block->hasInterpolationFlagCF(D3Q27System::S))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1+1, fbx2, fbx3, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1+1, fbx2, fbx3+0, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1+1, fbx2, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BSE);
      }
      if (block->hasInterpolationFlagCF(D3Q27System::BNW)&& !block->hasInterpolationFlagCF(D3Q27System::BN)&& !block->hasInterpolationFlagCF(D3Q27System::BW)&& !block->hasInterpolationFlagCF(D3Q27System::NW)&& !block->hasInterpolationFlagCF(D3Q27System::B)&& !block->hasInterpolationFlagCF(D3Q27System::N) && !block->hasInterpolationFlagCF(D3Q27System::W))
      {
         SPtr<Block3D> fblock00 = grid->getBlock(fbx1, fbx2+1, fbx3+0, level);
         SPtr<Block3D> fblock10;// = grid->getBlock(fbx1, fbx2+1, fbx3, level);
         SPtr<Block3D> fblock01;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);
         SPtr<Block3D> fblock11;// = grid->getBlock(fbx1, fbx2+1, fbx3+1, level);

         setInterpolationConnectors(fblock00, fblock10, fblock01, fblock11, block, D3Q27System::BNW);
      }

   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setInterpolationConnectors() - end");
}
//////////////////////////////////////////////////////////////////////////

void ConnectorBlockVisitor::setInterpolationConnectors(SPtr<Block3D> fblock00, SPtr<Block3D> fblock10, SPtr<Block3D> fblock01, SPtr<Block3D> fblock11, SPtr<Block3D> cBlock, int dir)
{
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setInterpolationConnectors(...) - start");
   int fblock00Rank = -999, fblock10Rank = -999, fblock01Rank = -999, fblock11Rank = -999;
   if (fblock00) fblock00Rank = fblock00->getRank();
   if (fblock01) fblock01Rank = fblock01->getRank();
   if (fblock10) fblock10Rank = fblock10->getRank();
   if (fblock11) fblock11Rank = fblock11->getRank();
   int cBlockRank = cBlock->getRank();

   LBMReal omegaF;
   if (fblock00) omegaF = LBMSystem::calcCollisionFactor(nu, fblock00->getLevel());
   if (fblock01) omegaF = LBMSystem::calcCollisionFactor(nu, fblock01->getLevel());
   if (fblock10) omegaF = LBMSystem::calcCollisionFactor(nu, fblock10->getLevel());
   if (fblock11) omegaF = LBMSystem::calcCollisionFactor(nu, fblock11->getLevel());
   LBMReal omegaC = LBMSystem::calcCollisionFactor(nu, cBlock->getLevel());
   iProcessor->setOmegas(omegaC, omegaF);

   InterpolationProcessorPtr cIProcessor(iProcessor->clone());
   InterpolationProcessorPtr fIProcessor00(iProcessor->clone());
   InterpolationProcessorPtr fIProcessor10(iProcessor->clone());
   InterpolationProcessorPtr fIProcessor01(iProcessor->clone());
   InterpolationProcessorPtr fIProcessor11(iProcessor->clone());

   CreateTransmittersHelper::TransmitterPtr senderCF00, receiverCF00,
                                             senderCF01, receiverCF01,
                                             senderCF10, receiverCF10,
                                             senderCF11, receiverCF11,
                                             senderFC00, receiverFC00,
                                             senderFC01, receiverFC01,
                                             senderFC10, receiverFC10,
                                             senderFC11, receiverFC11;

   if (fblock00) createTransmitters(cBlock, fblock00, dir, CreateTransmittersHelper::SW, senderCF00, receiverCF00, senderFC00, receiverFC00);
   if (fblock01) createTransmitters(cBlock, fblock01, dir, CreateTransmittersHelper::NW, senderCF01, receiverCF01, senderFC01, receiverFC01);
   if (fblock10) createTransmitters(cBlock, fblock10, dir, CreateTransmittersHelper::SE, senderCF10, receiverCF10, senderFC10, receiverFC10);
   if (fblock11) createTransmitters(cBlock, fblock11, dir, CreateTransmittersHelper::NE, senderCF11, receiverCF11, senderFC11, receiverFC11);

   if (cBlockRank == gridRank)
   {
      SPtr<Block3DConnector> connector = cFactory->createCoarseToFineConnector(cBlock,
         senderCF00, receiverCF00, senderCF01, receiverCF01,
         senderCF10, receiverCF10, senderCF11, receiverCF11,
         dir, cIProcessor);
      cBlock->setConnector(connector);
   }
   if (fblock00 && fblock00Rank == gridRank)
   {
      SPtr<Block3DConnector> connector = cFactory->createFineToCoarseConnector(fblock00,
         senderFC00, receiverFC00, dir, fIProcessor00, FineToCoarseBlock3DConnector::Type00);
      fblock00->setConnector(connector);
   }
   if (fblock01 && fblock01Rank == gridRank)
   {
      SPtr<Block3DConnector> connector = cFactory->createFineToCoarseConnector(fblock01,
         senderFC01, receiverFC01, dir, fIProcessor01, FineToCoarseBlock3DConnector::Type01);
      fblock01->setConnector(connector);
   }
   if (fblock10 && fblock10Rank == gridRank)
   {
      SPtr<Block3DConnector> connector = cFactory->createFineToCoarseConnector(fblock10,
         senderFC10, receiverFC10, dir, fIProcessor10, FineToCoarseBlock3DConnector::Type10);
      fblock10->setConnector(connector);
   }
   if (fblock11 && fblock11Rank == gridRank)
   {
      SPtr<Block3DConnector> connector = cFactory->createFineToCoarseConnector(fblock11,
         senderFC11, receiverFC11, dir, fIProcessor11, FineToCoarseBlock3DConnector::Type11);
      fblock11->setConnector(connector);
   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::setInterpolationConnectors(...) - end");
}
//////////////////////////////////////////////////////////////////////////
void ConnectorBlockVisitor::createTransmitters(SPtr<Block3D> cBlock, SPtr<Block3D> fBlock, int dir,
   CreateTransmittersHelper::IBlock ib,
   CreateTransmittersHelper::TransmitterPtr& senderCF,
   CreateTransmittersHelper::TransmitterPtr& receiverCF,
   CreateTransmittersHelper::TransmitterPtr& senderFC,
   CreateTransmittersHelper::TransmitterPtr& receiverFC)
{
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::createTransmitters(...) - start");
   CreateTransmittersHelper helper;
//   bool MPIpool = true;
//   bool orthogonal = false;
   int fBlockRank = fBlock->getRank();
   int cBlockRank = cBlock->getRank();
   if (fBlockRank == cBlockRank && fBlockRank == gridRank)
   {
      senderCF = receiverFC = CreateTransmittersHelper::TransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
      senderFC = receiverCF = CreateTransmittersHelper::TransmitterPtr(new TbLocalTransmitter< CbVector< LBMReal > >());
   }
   else if (cBlockRank == gridRank)
   {
      helper.createTransmitters(cBlock, fBlock, dir, ib, senderCF, receiverCF, comm, CreateTransmittersHelper::MPI);
   }
   else if (fBlockRank == gridRank)
   {
      helper.createTransmitters(fBlock, cBlock, dir, ib, senderFC, receiverFC, comm, CreateTransmittersHelper::MPI);
   }
   UBLOG(logDEBUG5, "ConnectorBlockVisitor::createTransmitters(...) - end");
}

