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

      if (comm->isRoot())
         UBLOG(logINFO, "RestartDemObjectsCoProcessor::write step: " << istep);

      write(istep);
   }
}

void RestartDemObjectsCoProcessor::restart(double step)
{
   if (comm->isRoot())
      UBLOG(logINFO, "RestartDemObjectsCoProcessor::read step: " << (int)step);

   read((int)step);
}

void RestartDemObjectsCoProcessor::write(int step)
{
   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor::write start ");
   std::vector<double> p;

   demCoProcessor->getObjectsPropertiesVector(p);

   //TODO implement getherv 
   std::vector<double> rvalues;
   comm->allGather(p, rvalues);

   if (comm->isRoot())
   {
      std::map< int, std::vector< double> > infMap;
      int size =  (int)rvalues.size();
      for (int i = 0; i < size; i += 7)
      {
         std::vector< double> infVector(6);
         for (int j = 0; j < 6; j ++)
         {
            infVector[j] = rvalues[i+1+j];
         }
         infMap.insert(std::make_pair((int)rvalues[i], infVector));
      }
      std::vector< double> wvalues;
      typedef std::map< int, std::vector< double> >::iterator it_type;
      for (it_type iterator = infMap.begin(); iterator != infMap.end(); iterator++) 
      {
         // iterator->first = key
         // iterator->second = value
         std::vector<double>::iterator it = wvalues.end();
         it = wvalues.insert(it, iterator->second.begin(), iterator->second.end());
      }
      std::string subfolder = "dem_cp_"+UbSystem::toString(step);
      std::string filePath =  path+"/dem_cp/"+subfolder+"/dem_cp.bin";
      UbFileOutputBinary fo(filePath);
      fo.writeInteger((int)wvalues.size());
      fo.writeVector<double>(wvalues);
      UBLOG(logINFO, "RestartDemObjectsCoProcessor::write number of objects = " << wvalues.size()/6);
   }
   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor::write stop ");
}

void RestartDemObjectsCoProcessor::read(int step)
{
   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor::read start ");
   std::vector<double> p;

   if (comm->isRoot())
   {
      std::string subfolder = "dem_cp_"+UbSystem::toString(step);
      std::string filePath =  path+"/dem_cp/"+subfolder+"/dem_cp.bin";
      UbFileInputBinary fi(filePath);
      int size = fi.readInteger();
      p.resize(size);
      fi.readVector<double>(p);
   }
   comm->broadcast(p);

   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor::read number of objects = " << p.size()/6);

   createDemObjectsCoProcessor->clearGeoObjects();

   int size =  (int)p.size();

   for (int i = 0; i < size; i += 6)
   {
      SPtr<GbObject3D> sphere(new GbSphere3D(p[i], p[i+1], p[i+2], radius));
      createDemObjectsCoProcessor->addGeoObject(sphere, Vector3D(p[i+3], p[i+4], p[i+5]));
   }

   createDemObjectsCoProcessor->createGeoObjects();

   createDemObjectsCoProcessor->clearGeoObjects();

   if (comm->isRoot()) UBLOG(logINFO, "RestartDemObjectsCoProcessor::read stop ");
}
