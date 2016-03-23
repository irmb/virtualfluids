#include "EmergencyExitPostprocessor.h"
#include <basics/utilities/UbFileOutputASCII.h>
#include <basics/utilities/UbFileInputASCII.h>

EmergencyExitPostprocessor::EmergencyExitPostprocessor( Grid3DPtr grid, UbSchedulerPtr s, 
                                                        const std::string& path, 
                                                        RestartPostprocessorPtr rp, CommunicatorPtr comm ) :
                                                        Postprocessor(grid, s),
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
EmergencyExitPostprocessor::~EmergencyExitPostprocessor()
{

}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitPostprocessor::update( double step )
{
   if(scheduler->isDue(step) )
      collectPostprocessData(step);

   UBLOG(logDEBUG3, "EmergencyExitPostprocessor::update:" << step);
}

void EmergencyExitPostprocessor::collectPostprocessData( double step )
{
   if(readMetafile())
   {
      rp->doCheckPoint((int)step);
      if(comm->getProcessID() == comm->getRoot()) UBLOG(logINFO,"EmergencyExitPostprocessor save step: " << step);
      comm->barrier();
      exit(EXIT_SUCCESS);
   }
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitPostprocessor::writeMetafile(int status )
{
   UbFileOutputASCII out(metafile);
   out.writeBool(false);
}
//////////////////////////////////////////////////////////////////////////
bool EmergencyExitPostprocessor::readMetafile()
{
   UbFileInputASCII in(metafile);
   return in.readBool();
}
//////////////////////////////////////////////////////////////////////////
void EmergencyExitPostprocessor::checkMetafile()
{
   std::ifstream file(metafile.c_str()); 
   if (!file.is_open()) 
   {
      writeMetafile(false);
      return;
   }
   file.close();
}
                                                       

