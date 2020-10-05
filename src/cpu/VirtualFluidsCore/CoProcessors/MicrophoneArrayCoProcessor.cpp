#include "MicrophoneArrayCoProcessor.h"
#include "Vector3D.h"
#include "Grid3D.h"
#include "ILBMKernel.h"
#include "Communicator.h"
#include "Block3D.h"
#include "DistributionArray3D.h"
#include "DataSet3D.h"
#include "D3Q27System.h"
#include "UbScheduler.h"
#include "BCProcessor.h"
#include "BCArray3D.h"
#include <sstream>

MicrophoneArrayCoProcessor::MicrophoneArrayCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string & path, SPtr<Communicator> comm) : CoProcessor(grid, s), path(path), comm(comm)
{
   count = 0;
   micID = 0;
}

MicrophoneArrayCoProcessor::~MicrophoneArrayCoProcessor()
{
}

void MicrophoneArrayCoProcessor::process(double step)
{
   if (microphones.size() > 0)
   {
      collectData(step);

      if (scheduler->isDue(step))
         writeFile(step);
   }

   UBLOG(logDEBUG3, "MicrophoneArrayCoProcessor::process:" << step);
}

bool MicrophoneArrayCoProcessor::addMicrophone(Vector3D coords)
{
   micID++;
//   UbTupleInt3 blockIndexes = grid->getBlockIndexes(coords[0], coords[1], coords[2]);

   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   for (int level = minInitLevel; level<=maxInitLevel; level++)
   {
      UbTupleInt3 blockIndexes = grid->getBlockIndexes(coords[0], coords[1], coords[2], level);
      SPtr<Block3D> block = grid->getBlock(val<1>(blockIndexes), val<2>(blockIndexes), val<3>(blockIndexes), level);
      if (block)
      {
         SPtr<ILBMKernel> kernel = block->getKernel();
         if (kernel)
         {
            SPtr<BCArray3D> bcarray = kernel->getBCProcessor()->getBCArray();
            UbTupleInt3 nodes = grid->getNodeIndexes(block, coords[0], coords[1], coords[2]);
            if (!bcarray->isUndefined(val<1>(nodes), val<2>(nodes), val<3>(nodes)))
            {

               if (kernel->getCompressible())
               {
                  calcMacros = &D3Q27System::calcCompMacroscopicValues;
               }
               else
               {
                  calcMacros = &D3Q27System::calcIncompMacroscopicValues;
               }
               SPtr<Mic> mic(new Mic);
               mic->id = micID;
               mic->distridution = kernel->getDataSet()->getFdistributions();
               mic->nodeIndexes = grid->getNodeIndexes(block, coords[0], coords[1], coords[2]);
               microphones.push_back(mic);

               strVector.push_back(SPtr<std::stringstream>(new std::stringstream));

               std::string fname = path+"/mic/mic_"+UbSystem::toString(micID)+".csv";
               std::ofstream ostr;
               ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
               if (!ostr)
               {
                  ostr.clear();
                  std::string path = UbSystem::getPathFromString(fname);
                  if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
                  if (!ostr) throw UbException(UB_EXARGS, "couldn't open file "+fname);
               }
               ostr << "#microphone position: " << coords[0] << "; " << coords[1] << "; " << coords[2] << "; " << "\n";
               ostr.close();
               return true;
            }
         }
      }
   }
   return false;
}

void MicrophoneArrayCoProcessor::collectData(double step)
{
   for (int i = 0; i < microphones.size(); i++ )
   {
      LBMReal f[D3Q27System::ENDF+1];
      microphones[i]->distridution->getDistribution(f, val<1>(microphones[i]->nodeIndexes), val<2>(microphones[i]->nodeIndexes), val<3>(microphones[i]->nodeIndexes));
      LBMReal vx1, vx2, vx3, rho;
      calcMacros(f, rho, vx1, vx2, vx3);
      *strVector[i] << step << ';' << rho << '\n';
   }
}

void MicrophoneArrayCoProcessor::writeFile(double step)
{
   for (int i = 0; i < microphones.size(); i++)
   {
      std::string fname = path+"/mic/mic_"+UbSystem::toString(microphones[i]->id)+".csv";
      std::ofstream ostr;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if (!ostr)
      {
         ostr.clear();
         std::string path = UbSystem::getPathFromString(fname);
         if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file "+fname);
      }
      ostr << strVector[i]->str();
      ostr.close();
      strVector[i] = SPtr<std::stringstream>(new std::stringstream);
   }
}