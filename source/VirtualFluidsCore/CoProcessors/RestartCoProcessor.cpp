#include "RestartCoProcessor.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbFileInputASCII.h>

#include "BoostSerializationClassExportHelper.h"

#include <MemoryUtil.h>

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

      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"RestartPostprocessor save step: " << step);

      UBLOG(logDEBUG3, "RestartPostprocessor::update:" << step);
   }   
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::addCoProcessor( PostprocessorPtr p )
{
   postprocessors.push_back(p);
}
//////////////////////////////////////////////////////////////////////////
PostprocessorPtr RestartCoProcessor::getCoProcessor(int index)
{
   return postprocessors[index];
}
//////////////////////////////////////////////////////////////////////////
std::vector<PostprocessorPtr> RestartCoProcessor::getPostprocessors()
{
   return postprocessors;
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::doCheckPoint(int step)
{
   UBLOG(logDEBUG3,"Save check point - start");

   int pid = comm->getProcessID();
   std::string filename = path + "/checkpoint" + UbSystem::toString(step) + "/checkpoint" + UbSystem::toString(pid) + "_" + UbSystem::toString(step);

   if (archiveType == TXT)
   {
      saveTxtArchive(filename + ".txt");
   } 
   else if(archiveType == BINARY)
   {
      saveBinArchive(filename + ".bin");
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

      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"RestartPostprocessor restart step: " << restartStep);

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
void RestartCoProcessor::saveTxtArchive(std::string filename)
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

   int psize = (int)postprocessors.size(); 
   oa << psize;

   oa.register_type<PostprocessorPtr>();
   BOOST_FOREACH(PostprocessorPtr pp, postprocessors)
   {
      oa << pp;
   }
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

   ia.register_type<PostprocessorPtr>();
   for (int i = 0; i < psize; i++)
   {
      PostprocessorPtr pp;
      ia >> pp;
      pp->reconnect(grid);
      postprocessors.push_back(pp);
   }
}
//////////////////////////////////////////////////////////////////////////
void RestartCoProcessor::saveBinArchive( std::string filename )
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

   int psize = (int)postprocessors.size(); 
   oa << psize;

   oa.register_type<PostprocessorPtr>();
   BOOST_FOREACH(PostprocessorPtr pp, postprocessors)
   {
      oa << pp;
   }
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

   ia.register_type<PostprocessorPtr>();
   for (int i = 0; i < psize; i++)
   {
      PostprocessorPtr pp;
      ia >> pp;
      pp->reconnect(grid);
      postprocessors.push_back(pp);
   }
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



