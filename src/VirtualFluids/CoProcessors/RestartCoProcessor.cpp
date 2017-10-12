#include "RestartCoProcessor.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbFileInputASCII.h>

#include "MetisPartitioningGridVisitor.h"

#include "BoostSerializationClassExportHelper.h"

#include <MemoryUtil.h>
#include <vector>
#include <boost/foreach.hpp>

RestartCoProcessor::RestartCoProcessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, ArchiveType type) :
   CoProcessor(grid, s),
   grid(grid),
   path(path),
   archiveType(type),
   comm(comm)
{
   restartStep = 0;
   this->path = path + "/checkpoints";
   metafile = this->path + "/LastCP.txt";
   if (comm->getProcessID() == comm->getRoot())
   {
      checkMetafile();
   }
   comm->barrier();
   grid = restart();
}
//////////////////////////////////////////////////////////////////////////
RestartCoProcessor::RestartCoProcessor(Grid3DPtr& grid, UbSchedulerPtr s, CommunicatorPtr comm, const std::string& path, int restartStep, ArchiveType type) :
CoProcessor(grid, s),
grid(grid),
path(path),
archiveType(type),
comm(comm),
restartStep(restartStep)
{
   this->path = path + "/checkpoints";
   metafile = this->path + "/LastCP.txt";
   if (comm->getProcessID() == comm->getRoot())
   {
      writeMetafile(restartStep);
   }
   comm->barrier();
   grid = restart();
}
//////////////////////////////////////////////////////////////////////////
RestartCoProcessor::~RestartCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::process(double step)
{
   if(scheduler->isDue(step) && step != restartStep)
   {
      doCheckPoint(int(step));

      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"RestartCoProcessor save step: " << step);

      UBLOG(logDEBUG3, "RestartCoProcessor::update:" << step);
   }   
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::addCoProcessor( CoProcessorPtr p )
{
   CoProcessors.push_back(p);
}
//////////////////////////////////////////////////////////////////////////
CoProcessorPtr RestartCoProcessor::getCoProcessor(int index)
{
   return CoProcessors[index];
}
//////////////////////////////////////////////////////////////////////////
std::vector<CoProcessorPtr> RestartCoProcessor::getCoProcessors()
{
   return CoProcessors;
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::doCheckPoint(int step)
{
   UBLOG(logDEBUG3,"Save check point - start");

   int pid = comm->getProcessID();
   std::string filename = path + "/checkpoint" + UbSystem::toString(step) + "/checkpoint" + UbSystem::toString(pid) + "_" + UbSystem::toString(step);

   if (archiveType == TXT)
   {
      saveTxtArchive(filename + ".txt", grid);
   } 
   else if(archiveType == BINARY)
   {
      saveBinArchive(filename + ".bin", grid);
   }

   comm->barrier();

   if (comm->getProcessID() == comm->getRoot())
   {
      writeMetafile(step);
   }

   UBLOG(logDEBUG3,"Save check point - end");
}
//////////////////////////////////////////////////////////////////////////
Grid3DPtr RestartCoProcessor::restart()
{
   restartStep = readMetafile();

   comm->barrier();

   if (restartStep > 0)
   {
      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"Load check point - start");
      int pid = comm->getProcessID();
      std::string filename = path + "/checkpoint" + UbSystem::toString(restartStep) + "/checkpoint" + UbSystem::toString(pid) + "_" + UbSystem::toString(restartStep);

      if (archiveType == TXT)
      {
         loadTxtArchive(filename + ".txt");
      } 
      else if(archiveType == BINARY)
      {
         loadBinArchive(filename + ".bin");
      }

      this->reconnect(grid);
      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"Load check point - end");

      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"RestartCoProcessor restart step: " << restartStep);

      return grid;
   } 
   else
   {
      return grid;
   }

}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::acceptGridVisitors()
{
   BOOST_FOREACH(Grid3DVisitorPtr v, gridVisitors)
   {
      grid->accept(*(v.get()));
   }
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::acceptBlockVisitors()
{
   BOOST_FOREACH(Block3DVisitorPtr v, blockVisitors)
   {
      grid->accept(*(v.get()));
   }
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::addGridVisitor( Grid3DVisitorPtr v )
{
   gridVisitors.push_back(v);
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::addBlockVisitor( Block3DVisitorPtr v )
{
   blockVisitors.push_back(v);
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::saveTxtArchive( std::string filename, Grid3DPtr grid )
{
   std::ofstream file(filename.c_str()); 
   if(!file)
   { 
      file.clear(); 
      std::string path = UbSystem::getPathFromString(filename);
      if(path.size()>0){ UbSystem::makeDirectory(path); file.open(filename.c_str());}
      if(!file) throw UbException(UB_EXARGS,"couldn't open file "+filename);
   }
   boost::archive::text_oarchive oa(file);
   oa.register_type<Grid3D>();
   oa << grid;

   int psize = (int)CoProcessors.size(); 
   oa << psize;

   oa.register_type<CoProcessorPtr>();
   BOOST_FOREACH(CoProcessorPtr pp, CoProcessors)
   {
      oa << pp;
   }
   file.close();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::loadTxtArchive( std::string filename )
{
   std::ifstream file(filename.c_str()); 
   if (!file.is_open()) UB_THROW( UbException(UB_EXARGS,"Can not open check point file \"" + filename + "\""));
   boost::archive::text_iarchive ia(file);
   ia.register_type<Grid3D>();
   ia >> grid; 

   int psize;
   ia >> psize;

   ia.register_type<CoProcessorPtr>();
   for (int i = 0; i < psize; i++)
   {
      CoProcessorPtr pp;
      ia >> pp;
      pp->reconnect(grid);
      CoProcessors.push_back(pp);
   }
   file.close();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::saveBinArchive( std::string filename, Grid3DPtr grid )
{
   //impotent for binary archive add std::ios::binary
   std::ofstream file(filename.c_str(), std::ios::binary); 
   if(!file)
   { 
      file.clear(); 
      std::string path = UbSystem::getPathFromString(filename);
      if(path.size()>0){ UbSystem::makeDirectory(path); file.open(filename.c_str());}
      if(!file) throw UbException(UB_EXARGS,"couldn't open file "+filename);
   }
   boost::archive::binary_oarchive oa(file);
   oa.register_type<Grid3D>();
   oa << grid;

   int psize = (int)CoProcessors.size(); 
   oa << psize;

   oa.register_type<CoProcessorPtr>();
   BOOST_FOREACH(CoProcessorPtr pp, CoProcessors)
   {
      oa << pp;
   }
   file.close();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::loadBinArchive( std::string filename )
{
   //impotent for binary archive add std::ios::binary
   std::ifstream file(filename.c_str(), std::ios::binary); 
   if (!file.is_open()) UB_THROW( UbException(UB_EXARGS,"Can not open check point file \"" + filename + "\""));
   boost::archive::binary_iarchive ia(file);
   ia.register_type<Grid3D>();
   ia >> grid;

   int psize;
   ia >> psize;

   ia.register_type<CoProcessorPtr>();
   for (int i = 0; i < psize; i++)
   {
      CoProcessorPtr pp;
      ia >> pp;
      pp->reconnect(grid);
      CoProcessors.push_back(pp);
   }
   file.close();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::writeMetafile(int step )
{
   UbFileOutputASCII out(metafile);
   out.writeInteger(step);
}
//////////////////////////////////////////////////////////////////////////
int RestartCoProcessor::readMetafile()
{
   UbFileInputASCII in(metafile);
   return in.readInteger();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::checkMetafile()
{
   std::ifstream file(metafile.c_str()); 
   if (!file.is_open()) 
   {
      writeMetafile(0);
      return;
   }
   file.close();
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::writeDistributedGrid(Grid3DPtr sgrid, int numberOfProcesses)
{
   using namespace std;

   Grid3DVisitorPtr metisVisitor(new MetisPartitioningGridVisitor(comm, MetisPartitioningGridVisitor::LevelBased, D3Q27System::BSW, MetisPartitioner::RECURSIVE));
   boost::dynamic_pointer_cast<MetisPartitioningGridVisitor>(metisVisitor)->setNumberOfProcesses(numberOfProcesses);
   sgrid->accept(metisVisitor);

   int minInitLevel = sgrid->getCoarsestInitializedLevel();
   int maxInitLevel = sgrid->getFinestInitializedLevel();

   for (int i = 0; i<numberOfProcesses; i++)
   {
      UBLOG(logINFO, "Create dump file for rank " << i <<" - start");
      Grid3DPtr newGrid(new Grid3D());
      newGrid->setRank(i);
      newGrid->setDeltaX(sgrid->getDeltaX(0));
      newGrid->setNX1(sgrid->getNX1());
      newGrid->setNX2(sgrid->getNX2());
      newGrid->setNX3(sgrid->getNX3());
      newGrid->setCoordinateTransformator(sgrid->getCoordinateTransformator());
      UbTupleInt3 blockNX = sgrid->getBlockNX();
      newGrid->setBlockNX(val<1>(blockNX), val<2>(blockNX), val<3>(blockNX));
      Grid3D::Interactor3DSet interactors = sgrid->getInteractors();
      for(int inter=0; inter < interactors.size(); inter++)
      {
         newGrid->addInteractor(interactors[inter]);
      }

      for (int level = minInitLevel; level<=maxInitLevel; level++)
      {
         vector<Block3DPtr> blockVector;
         grid->getBlocks(level, blockVector);
         BOOST_FOREACH(Block3DPtr block, blockVector)
         {
            if (block)
            {
               if (block->getRank() == i)
               {
                  newGrid->addBlock(block);
               } 
               else
               {
                  Block3DPtr newBlock(new Block3D(block->getX1(), block->getX2(), block->getX3(), block->getLevel()));
                  newBlock->setRank(block->getRank());
                  newGrid->addBlock(newBlock);
               }
            }
         }
      }

      std::string filename = path+"/checkpoint"+UbSystem::toString(1)+"/checkpoint"+UbSystem::toString(i)+"_"+UbSystem::toString(1);

      if (archiveType==TXT)
      {
         saveTxtArchive(filename+".txt", newGrid);
      }
      else if (archiveType==BINARY)
      {
         saveBinArchive(filename+".bin", newGrid);
      }

      UBLOG(logINFO, "Create dump file for rank " << i <<" - end");
   }
   writeMetafile(1);
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::setArchiveType(ArchiveType type)
{
   archiveType = type;
}
//////////////////////////////////////////////////////////////////////////
RestartCoProcessor::ArchiveType RestartCoProcessor::getArchiveType()
{
   return archiveType;
}


