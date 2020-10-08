#include "EmergencyExitCoProcessor.h"
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbFileInputASCII.h>
#include "UbLogger.h"
#include "UbScheduler.h"
#include "Communicator.h"
#include "MPIIORestartCoProcessor.h"
#include "Grid3D.h"

EmergencyExitCoProcessor::EmergencyExitCoProcessor( SPtr<Grid3D> grid, SPtr<UbScheduler> s, 
                                                        const std::string& path, 
                                                        SPtr<MPIIORestartCoProcessor> rp, SPtr<Communicator> comm) :
                                                        CoProcessor(grid, s),
                                                        path(path),
                                                        rp(rp),
                                                        comm(comm)
{
   this->path = path + "/exit";
   metafile = this->path + "/stop.txt";
   if (comm->getProcessID() == comm->getRoot())
   {
      //checkMetafile();
      writeMetafile(false);
   }
   comm->barrier();
}
//////////////////////////////////////////////////////////////////////////
EmergencyExitCoProcessor::~EmergencyExitCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitCoProcessor::process( double step )
{
   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "EmergencyExitCoProcessor::update:" << step);
}

void EmergencyExitCoProcessor::collectData( double step )
{
   if(readMetafile())
   {
      rp->process((int)step);
      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"EmergencyExitCoProcessor save step: " << step);
      comm->barrier();
      exit(EXIT_SUCCESS);
   }
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitCoProcessor::writeMetafile(int  /*status*/ )
{
   UbFileOutputASCII out(metafile);
   out.writeBool(false);
}
//////////////////////////////////////////////////////////////////////////
bool EmergencyExitCoProcessor::readMetafile()
{
   UbFileInputASCII in(metafile);
   return in.readBool();
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitCoProcessor::checkMetafile()
{
   std::ifstream file(metafile.c_str()); 
   if (!file.is_open()) 
   {
      writeMetafile(false);
      return;
   }
   file.close();
}
                                                       

