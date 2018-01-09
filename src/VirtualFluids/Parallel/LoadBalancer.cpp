#if defined VF_MPI

#include "LoadBalancer.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>
#include "MPICommunicator.h"
#include "MathUtil.hpp"

//#include "BoostSerializationClassExportHelper.h"

LoadBalancer::LoadBalancer(Grid3DPtr grid, CommunicatorPtr comm, const int& endDir) : 
grid(grid),   
comm(comm),
endDir(endDir)
{
   mpi_comm = *((MPI_Comm*) comm->getNativeCommunicator());
   processID = comm->getProcessID();
   int size = comm->getNumberOfProcesses();
   converged.resize(size,0);
}
//////////////////////////////////////////////////////////////////////////
LoadBalancer::~LoadBalancer()
{

}
//////////////////////////////////////////////////////////////////////////
bool LoadBalancer::balance()
{
   if (!isConverged())
   {
      //std::cout << "is not converged" << std::endl;
      if(processID == 0)
         UBLOG(logINFO,"is not converged");
      collectData();
      collectNeighboursLoad();
      prepareToSendLoad();
      prepareToRecieveLoad();
      sendRecieveLoad();
      saveRecievedLoad();
      //MPI_Barrier(mpi_comm);
      //std::cout<< "MPI_Barrier: PID = " << processID <<std::endl;
      return false;
   }
   else
   {
      //std::cout << "is converged" << std::endl;
      if(processID == 0)
         UBLOG(logINFO,"is converged");
      return true;
   }
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::prepareToSendLoad()
{
   sendBuffer.clear();
   removeBlocksSend.resize(0);
   int conv = 1;
   for(NeighbourProcess neighbProzess : neighbourProcesses)
   {
      if (this->load != neighbProzess.load)
      {
         if(neighbProzess.load < this->load)
         {
            //std::cout<<"rank = " << comm->getProcessID() << " : send Abweichung =  " << (this->load/neighbProzess.load - 1.0)*100.0 <<std::endl;
            if ((static_cast<double>(this->load)/static_cast<double>(neighbProzess.load) - 1.0)*100.0 > 3.0)
            {
               //if(processID == 0)
                  UBLOG(logINFO,"rank = " << comm->getProcessID() << " : send Abweichung =  " << (this->load/neighbProzess.load - 1.0)*100.0);
               double  transLoad = (this->load - neighbProzess.load)/2;
               prepareToSendLoad(neighbProzess.processID, transLoad, conv);
            }
         }
      }
   }
   setConverged(conv);
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::prepareToSendLoad(const int& nProcessID, const double& transLoad, int& conv)
{
   if((itnMap = neighboursMap.find(nProcessID)) != neighboursMap.end())
   {
      TransferBuffer sbuf;
      std::vector<Block3DPtr> sendBlocks;
      BBlocks &boundaryBlocks = itnMap->second;
      TBlocks tBlocks;
      std::stringstream  ss_out;

      std::string str_out;

      for(LoadBalancer::BBlocks::value_type b : boundaryBlocks)
      {
         TransferBlock &tBlock = b.second;
         tBlocks.insert(std::make_pair(tBlock.weight, tBlock));
      }

      double lsum = 0;
      for(LoadBalancer::TBlocks::value_type t : tBlocks)
      {
         TransferBlock &tBlock = t.second;
         lsum += tBlock.load;

         if(lsum >= transLoad)
            break;

         sendBlocks.push_back(tBlock.block);
         for(BBlocksMap::value_type &bbm : neighboursMap)
            bbm.second.erase(tBlock.block);

         //if(lsum >= transLoad)
         //   break;
      }
      if (sendBlocks.size() > 0)
      {
         boost::archive::text_oarchive oa(ss_out);
         oa.register_type<Block3D>();

         for(Block3DPtr b : sendBlocks)
         {
            oa << b;
         }

         sbuf.sbuffer = ss_out.str();
         sbuf.sbufferSize = static_cast<int> (sbuf.sbuffer.length()); 
         sbuf.blocksSize = static_cast<int>(sendBlocks.size());

         sbuf.comFlag = 1;
         sendBuffer.insert(std::make_pair(nProcessID, sbuf));

         for(Block3DPtr b : sendBlocks)
         {
            b->deleteKernel();
            b->setRank(nProcessID);
            b->deleteConnectors();
            //b->deleteInterpolationConnectors();

            removeBlocksSend.push_back(b->getX1());
            removeBlocksSend.push_back(b->getX2());
            removeBlocksSend.push_back(b->getX3());
            removeBlocksSend.push_back(b->getLevel());
            removeBlocksSend.push_back(b->getRank());

            //std::cout<< "Send: PID = " << processID << " block: " <<b->toString()<<std::endl;
         }
         conv *= 0;
      } 
      else
      {
         sbuf.comFlag = 0;
         sendBuffer.insert(std::make_pair(nProcessID, sbuf));
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::prepareToRecieveLoad()
{
   receiveBuffer.clear();
   int conv = 1;
   for(NeighbourProcess neighbProzess : neighbourProcesses)
   {
      if (this->load != neighbProzess.load)
      {
         if(neighbProzess.load > this->load)
         {
            //std::cout<<"rank = " << comm->getProcessID() << " : recieve Abweichung =  " << (static_cast<double>(neighbProzess.load)/static_cast<double>(this->load) - 1.0)*100.0 <<std::endl;
            if ((static_cast<double>(neighbProzess.load)/static_cast<double>(this->load) - 1.0)*100.0 > 3.0)
            {
               //if(processID == 0)
                  UBLOG(logINFO,"rank = " << comm->getProcessID() << " : recieve Abweichung =  " << (static_cast<double>(neighbProzess.load)/static_cast<double>(this->load) - 1.0)*100.0);
               conv *= 0;
               TransferBuffer rbuf;
               receiveBuffer.insert(std::make_pair(neighbProzess.processID, rbuf));
            }
         }
      }
   }
   setConverged(conv);
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::sendRecieveLoad()
{
   std::vector<MPI_Request> request;
   request.resize(0);
   int rcount = 0;
   for(TransferBufferMap::value_type &tb : receiveBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &rbuf = tb.second;
      request.push_back(0);
      MPI_Irecv(&rbuf.comFlag, 1, MPI_INT, nProcessID, 0, mpi_comm, &request[rcount]);
      rcount++;
   }
   for(TransferBufferMap::value_type &tb : sendBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &sbuf = tb.second;
      request.push_back(0);
      MPI_Isend(&sbuf.comFlag, 1, MPI_INT, nProcessID, 0, mpi_comm, &request[rcount]);
      rcount++;
   }
   if(request.size() > 0)
      MPI_Waitall(static_cast<int>(request.size()), &request[0], MPI_STATUSES_IGNORE);

   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "sendRecieveLoad comFlag" );
   //////////////////////////////////////////////////////////////////////////
   request.resize(0);
   rcount = 0;
   for(TransferBufferMap::value_type &tb : receiveBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &rbuf = tb.second;
      if (rbuf.comFlag == 1)
      {
         request.push_back(0);
         MPI_Irecv(&rbuf.blocksSize, 1, MPI_INT, nProcessID, 0, mpi_comm, &request[rcount]);
         rcount++;
      }
   }
   for(TransferBufferMap::value_type &tb : sendBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &sbuf = tb.second;
      if (sbuf.comFlag == 1)
      {
         request.push_back(0);
         MPI_Isend(&sbuf.blocksSize, 1, MPI_INT, nProcessID, 0, mpi_comm, &request[rcount]);
         rcount++;
      }
   }
   if(request.size() > 0)
      MPI_Waitall(static_cast<int>(request.size()), &request[0], MPI_STATUSES_IGNORE);
   
   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "sendRecieveLoad blocksSize" );
   //////////////////////////////////////////////////////////////////////////
   request.resize(0);
   rcount = 0;
   for(TransferBufferMap::value_type &tb : receiveBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &rbuf = tb.second;
      if (rbuf.comFlag == 1)
      {
         request.push_back(0);
         MPI_Irecv(&rbuf.sbufferSize, 1, MPI_INT, nProcessID, 0, mpi_comm, &request[rcount]);
         rcount++;
      }
   }
   for(TransferBufferMap::value_type &tb : sendBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &sbuf = tb.second;
      if (sbuf.comFlag == 1)
      {
         request.push_back(0);
         MPI_Isend(&sbuf.sbufferSize,1,MPI_INT,nProcessID,0,mpi_comm,&request[rcount]);
         rcount++;
      }
   }
   if(request.size() > 0)
      MPI_Waitall(static_cast<int>(request.size()), &request[0], MPI_STATUSES_IGNORE);

   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "sendRecieveLoad sbufferSize" );
   //////////////////////////////////////////////////////////////////////////
   request.resize(0);
   rcount = 0;
   for(TransferBufferMap::value_type &tb : receiveBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &rbuf = tb.second;
      if (rbuf.comFlag == 1)
      {
         request.push_back(0);
         rbuf.sbuffer.resize(rbuf.sbufferSize);
         MPI_Irecv((char *)rbuf.sbuffer.c_str(),rbuf.sbufferSize,MPI_CHAR,nProcessID,0,mpi_comm,&request[rcount]);
         rcount++;
      }
   }
   for(TransferBufferMap::value_type &tb : sendBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &sbuf = tb.second;
      if (sbuf.comFlag == 1)
      {
         request.push_back(0);
         MPI_Isend((char *)sbuf.sbuffer.c_str(),sbuf.sbufferSize,MPI_CHAR,nProcessID,0,mpi_comm,&request[rcount]);
         rcount++;
      }
   }
   if(request.size() > 0)
      MPI_Waitall(static_cast<int>(request.size()), &request[0], MPI_STATUSES_IGNORE);

   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "sendRecieveLoad sbuffer" );
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::sendRecieveChanges()
{
   std::vector<int> removeBlocksRecieve;

   //////////////////////////////////////////////////////////////////////////
   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "sendRecieveChanges" );
   
   comm->allGather(removeBlocksSend, removeBlocksRecieve);
   
   //////////////////////////////////////////////////////////////////////////
   //UBLOG(logINFO,"rank = " << comm->getProcessID() << "allGatherInts removeBlocksRecieve" );

   if(removeBlocksRecieve.size() >= 4)
   {
      for (int i = 0; i < removeBlocksRecieve.size(); i+=5)
      {
         Block3DPtr block = grid->getBlock(removeBlocksRecieve[i], removeBlocksRecieve[i+1], removeBlocksRecieve[i+2], removeBlocksRecieve[i+3]);
         block->setRank(removeBlocksRecieve[i+4]);
         //std::cout<< "RecieveChanges: PID = " << processID << " block: " <<block->toString()<<std::endl;
      }
   }
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::saveRecievedLoad()
{
   for(TransferBufferMap::value_type &tb : receiveBuffer)
   {
      int nProcessID = tb.first;
      TransferBuffer &rbuf = tb.second;

      if (rbuf.comFlag == 1)
      {
         std::stringstream  ss_in;
         ss_in.str(rbuf.sbuffer);

         boost::archive::text_iarchive ia(ss_in);
         ia.register_type<Block3D>();

         for(int i = 0; i < rbuf.blocksSize; i++)
         {
            Block3DPtr block;
            ia >> block;
            block->setRank(this->processID);
            grid->replaceBlock(block);
            //std::cout<< "Receive: PID = " << processID << " block: " <<block->toString()<<std::endl;
         }
      }
   }
   sendRecieveChanges();
}
//////////////////////////////////////////////////////////////////////////
bool LoadBalancer::isConverged()
{
   int conv = 1;
   for(int c : converged)
      conv *= c;
   
   if(conv)
      return true;
   else
      return false;
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::setConverged(int value)
{
   //UBLOG(logINFO,"rank = " << comm->getProcessID() << " setConverged start" );
   converged.resize(comm->getNumberOfProcesses());
   MPI_Allgather( &value, 1, MPI_INT, &converged[0], 1, MPI_INT, mpi_comm);
   //UBLOG(logINFO,"rank = " << comm->getProcessID() << " setConverged end" );
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::collectData()
{
   neighboursMap.clear();

   int weightScaleFactor = 8;

   int blockRank = 0;
   this->load = 0;

   int gridRank = grid->getRank();

   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      std::vector<Block3DPtr> blockVector;
      grid->getBlocks(level, gridRank, true, blockVector);
      for(Block3DPtr block : blockVector)
      {
         if (block)
         {
            this->load += block->getWorkLoad();
            blockRank  = block->getRank();
            
            //search  for parent block
            Block3DPtr parentBlock = grid->getSuperBlock(block);
            if(parentBlock)
            {
               int parentRank = parentBlock->getRank();
               if(parentRank != blockRank && parentBlock->isActive())
               {
                  block->addWeight(parentRank, weightScaleFactor);
                  addTransferBlock(block, parentRank);
               }
            }

            //search for child blocks
            std::vector<Block3DPtr> childBlocks;
            grid->getSubBlocks(block, 1, childBlocks);
            for(Block3DPtr b : childBlocks)
            {
               int childRank = b->getRank();
               if(childRank != blockRank && b->isActive())
               {
                  block->addWeight(childRank, weightScaleFactor);
                  addTransferBlock(block, childRank);
               }
            }

            //search for neighbor blocks
            for( int dir = 0; dir <= endDir; dir++)
            { 
               Block3DPtr neighBlock = (grid->getNeighborBlock(dir, block->getX1(), block->getX2(), block->getX3(), block->getLevel()));
               if(neighBlock)
               {
                  int neighBlockRank = neighBlock->getRank();
                  if(blockRank != neighBlockRank && neighBlock->isActive())
                  {
                     addTransferBlock(block, neighBlockRank);
                  }
               }
            }
         }
      }
   }

}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::collectNeighboursLoad()
{
   loadThreshold = 0;
   neighbourProcesses.resize(0);
   for(BBlocksMap::value_type b :  neighboursMap)
   {
      int neighbourProcessID = b.first;
      NeighbourProcess nProc;
      nProc.processID = neighbourProcessID;
      neighbourProcesses.push_back(nProc);
   }
   //std::cout<<"rank = " << comm->getProcessID() << " : this->load =  " << this->load <<std::endl;
   //if(processID == 0)
      UBLOG(logINFO, "rank = " << comm->getProcessID() << " : this->load =  " << this->load);

   std::vector<MPI_Request> request;
   request.resize(0);
   int rcount = 0;

   for(NeighbourProcess &np : neighbourProcesses)
   {
      request.push_back(0);
      MPI_Irecv(&np.load, 1, MPI_DOUBLE, np.processID, 0, mpi_comm, &request[rcount]);
      rcount++;
      //UBLOG(logINFO,"rank = " << comm->getProcessID() << " : MPI_Irecv from " << np.processID);
      request.push_back(0);
      MPI_Isend(&load, 1, MPI_DOUBLE, np.processID, 0, mpi_comm, &request[rcount]);
      rcount++;
      //UBLOG(logINFO,"rank = " << comm->getProcessID() << " : MPI_Isend to " << np.processID);
   }
   
   if(request.size() > 0)
      MPI_Waitall(static_cast<int>(request.size()), &request[0], MPI_STATUSES_IGNORE);

   for(NeighbourProcess &np : neighbourProcesses)
   {
      loadThreshold += np.load;
      //std::cout<<"rank = " << comm->getProcessID() << " : neighbor process load =  " << np.load <<std::endl;
      //UBLOG(logINFO,"rank = " << comm->getProcessID() << " : neighbor process load =  " << np.load);
   }
   loadThreshold += this->load;
   //int allLoad = loadThreshold;
   loadThreshold /= (static_cast<int>(neighbourProcesses.size()) + 1);
   //loadThresholdP = MathUtil::cint(static_cast<double>(loadThreshold)*100.0/static_cast<double>(allLoad));
   //loadP = MathUtil::cint(static_cast<double>(load)*100.0/static_cast<double>(allLoad));

   //std::cout<<"rank = " << comm->getProcessID() << " : this->loadThresholdP =  " << this->loadThresholdP <<std::endl;
   //std::cout<<"rank = " << comm->getProcessID() << " : this->loadP =  " << this->loadP <<std::endl;
   //std::cout<<"rank = " << comm->getProcessID() << " : allLoad =  " << allLoad <<std::endl;
}
//////////////////////////////////////////////////////////////////////////
void LoadBalancer::addTransferBlock(Block3DPtr block, int neighBlockRank)
{
   TransferBlock tBlock;
   tBlock.block = block;
   tBlock.load = block->getWorkLoad();
   tBlock.weight = block->getWeight(neighBlockRank);

   if((itnMap = neighboursMap.find(neighBlockRank)) == neighboursMap.end())
   {
      BBlocks boundaryBlocks;
      boundaryBlocks.insert(std::make_pair(block, tBlock));
      neighboursMap.insert(std::make_pair(neighBlockRank, boundaryBlocks));
   }
   else
   {
      BBlocks &boundaryBlocks = itnMap->second;
      boundaryBlocks.insert(std::make_pair(block, tBlock));
   }
}
//////////////////////////////////////////////////////////////////////////
#endif
