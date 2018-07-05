#include "RestartDemObjectsCoProcessor.h"

#include "Vector3D.h"
#include "Communicator.h"
#include "UbScheduler.h"
#include "Grid3D.h"
#include "UbSystem.h"
#include "GbSphere3D.h"
#include "DemCoProcessor.h"
#include "UbFileInputBinary.h"
#include "UbFileOutputBinary.h"
#include "CreateDemObjectsCoProcessor.h"

RestartDemObjectsCoProcessor::RestartDemObjectsCoProcessor()
{
}

RestartDemObjectsCoProcessor::RestartDemObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string & path, SPtr<DemCoProcessor> demCoProcessor, SPtr<CreateDemObjectsCoProcessor> createDemObjectsCoProcessor, double radius, SPtr<Communicator> comm)  : CoProcessor(grid, s), path(path), demCoProcessor(demCoProcessor), createDemObjectsCoProcessor(createDemObjectsCoProcessor), radius(radius), comm(comm)
{
}

void RestartDemObjectsCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
   {
      int istep = static_cast<int>(step);

      write(istep);

     if (comm->isRoot())
        UBLOG(logINFO, "RestartDemObjectsCoProcessor write step: " << istep);
   }
}

void RestartDemObjectsCoProcessor::restart(double step)
{
   read((int)step);
}

void RestartDemObjectsCoProcessor::write(int step)
{
   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor write:start ");
   std::vector<double> p;

   demCoProcessor->getObjectsPropertiesVector(p);

   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor size p: " << p.size());

   std::vector<double> rvalues;
   comm->allGather(p, rvalues);

   if (comm->isRoot())
   {
      std::string subfolder = "dem_cp_"+UbSystem::toString(step);
      std::string filePath =  path+"/dem_cp/"+subfolder+"/dem_cp.bin";
      UbFileOutputBinary fo(filePath);
      fo.writeInteger((int)rvalues.size());
      fo.writeVector<double>(rvalues);
      UBLOG(logINFO, "RestartDemObjectsCoProcessor size: " << rvalues.size());
   }
   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor write:stop ");
}

void RestartDemObjectsCoProcessor::read(int step)
{
   std::vector<double> p;

   if (comm->isRoot())
   {
      std::string subfolder = "dem_cp_"+UbSystem::toString(step);
      std::string filePath =  path+"/dem_cp/"+subfolder+"/dem_cp.bin";
      UbFileInputBinary fi(filePath);
      int size = fi.readInteger();
      p.resize(size);
      fi.readVector<double>(p);
      if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor read size p: " << p.size());
   }
   comm->broadcast(p);

   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor read size p broadcast: " << p.size());

   createDemObjectsCoProcessor->clearGeoObjects();

   int size =  (int)p.size();

   for (int i = 0; i < size; i += 6)
   {
      SPtr<GbObject3D> sphere(new GbSphere3D(p[i], p[i+1], p[i+2], radius));
      createDemObjectsCoProcessor->addGeoObject(sphere, Vector3D(p[i+3], p[i+4], p[i+5]));
   }

   createDemObjectsCoProcessor->createGeoObjects();

   createDemObjectsCoProcessor->clearGeoObjects();

}
