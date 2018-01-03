#include "WriteObjectsCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"

#include "basics/writer/WbWriterVtkXmlBinary.h"
#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "GbSphere3D.h"
#include "Communicator.h"
#include "UbScheduler.h"
#include "Grid3D.h"

WriteObjectsCoProcessor::WriteObjectsCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteObjectsCoProcessor::WriteObjectsCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
   const std::string& path, CommunicatorPtr comm)
   : CoProcessor(grid, s),
   path(path),
   comm(comm)
{
    this->path += "/geo/objects/sphere";
}


void WriteObjectsCoProcessor::addGbObject(GbSphere3DPtr sphere)
{
    objects.push_back(sphere);
}

void WriteObjectsCoProcessor::process(double step)
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
       std::string outFilename = WbWriterVtkXmlBinary::getInstance()->writeTriangles(path + std::to_string(stepInt), nodes, triangles);

   }
     

}
