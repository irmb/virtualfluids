#include "WriteDemObjectsCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"

#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "GbSphere3D.h"
#include "Communicator.h"
#include "UbScheduler.h"
#include "Grid3D.h"
#include "UbSystem.h"

WriteDemObjectsCoProcessor::WriteDemObjectsCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteDemObjectsCoProcessor::WriteDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, SPtr<DemCoProcessor> demCoProcessor, SPtr<Communicator> comm)
   : CoProcessor(grid, s),
   path(path),
   demCoProcessor(demCoProcessor),
   comm(comm)
{
    this->path += "/geo/objects/sphere";
}
//////////////////////////////////////////////////////////////////////////
void WriteDemObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
       std::vector<UbTupleFloat3> nodes;
       std::vector<UbTupleInt3>   triangles;

       for (size_t i = 0; i < objects.size(); i++)
       {
           objects[i]->addSurfaceTriangleSet(nodes, triangles);
  

       }
       int stepInt = (int)step;
       std::string outFilename = WbWriterVtkXmlBinary::getInstance()->writeTriangles(path + UbSystem::toString(stepInt), nodes, triangles);


       std::string pfilePath, partPath, subfolder, cfilePath;

       subfolder = "mq"+UbSystem::toString(istep);
       pfilePath =  path+"/mq/"+subfolder;
       cfilePath =  path+"/mq/mq_collection";
       partPath = pfilePath+"/mq"+UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep);


       std::string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);
       size_t found=partName.find_last_of("/");
       std::string piece = partName.substr(found+1);
       piece = subfolder + "/" + piece;

       std::vector<std::string> cellDataNames;
       std::vector<std::string> pieces = comm->gather(piece);
       if (comm->getProcessID() == comm->getRoot())
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
          UBLOG(logINFO, "WriteMacroscopicQuantitiesCoProcessor step: " << istep);
       }

   }
     

}
