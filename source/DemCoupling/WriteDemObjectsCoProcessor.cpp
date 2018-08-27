#include "WriteDemObjectsCoProcessor.h"

#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "Communicator.h"
#include "UbScheduler.h"
#include "Grid3D.h"
#include "UbSystem.h"
#include "DemCoProcessor.h"

WriteDemObjectsCoProcessor::WriteDemObjectsCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteDemObjectsCoProcessor::WriteDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<DemCoProcessor> demCoProcessor, SPtr<Communicator> comm)
   : CoProcessor(grid, s),
   path(path),
   writer(writer),
   demCoProcessor(demCoProcessor),
   comm(comm)
{
    
}
//////////////////////////////////////////////////////////////////////////
void WriteDemObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
       std::vector<UbTupleFloat3> nodes;
       std::vector<UbTupleInt3>   triangles;

       demCoProcessor->addSurfaceTriangleSet(nodes, triangles);

       int istep = static_cast<int>(step);

       std::string pfilePath, partPath, subfolder, cfilePath;

       subfolder = "dem"+UbSystem::toString(istep);
       pfilePath =  path+"/dem/"+subfolder;
       cfilePath =  path+"/dem/dem_collection";
       partPath = pfilePath+"/dem"+UbSystem::toString(comm->getProcessID())+ "_" + UbSystem::toString(istep);


       std::string partName = writer->writeTriangles(partPath, nodes, triangles);
       size_t found=partName.find_last_of("/");
       std::string piece = partName.substr(found+1);
       piece = subfolder + "/" + piece;

       std::vector<std::string> datanames;
       std::vector<std::string> cellDataNames;
       std::vector<std::string> pieces = comm->gather(piece);
       if (comm->isRoot())
       {
          std::string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
          found=pname.find_last_of("/");
          piece = pname.substr(found+1);

          std::vector<std::string> filenames;
          filenames.push_back(piece);
          if (step == CoProcessor::scheduler->getMinBegin())
          {
             WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
          }
          else
          {
             WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
          }
          UBLOG(logINFO, "WriteDemObjectsCoProcessor step: " << istep);
       }
   }
}
