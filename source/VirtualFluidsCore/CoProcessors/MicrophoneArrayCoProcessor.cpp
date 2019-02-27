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

MicrophoneArrayCoProcessor::MicrophoneArrayCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string & path, SPtr<Communicator> comm) : CoProcessor(grid, s), path(path), comm(comm)
{
   count = 0;
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

void MicrophoneArrayCoProcessor::addMicrophone(Vector3D coords)
{
   UbTupleInt3 blockIndexes = grid->getBlockIndexes(coords[0], coords[1], coords[2]);

   //int gridRank = comm->getProcessID();
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
               mic->distridution = kernel->getDataSet()->getFdistributions();
               mic->nodeIndexes = grid->getNodeIndexes(block, coords[0], coords[1], coords[2]);
               microphones.push_back(mic);
               values.resize((microphones.size()+1)*static_cast<int>(scheduler->getMinStep()));

               std::string fname = path+"/mic/mic_"+UbSystem::toString(comm->getProcessID())+".csv";
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
               return;
            }
         }
      }
   }
}

void MicrophoneArrayCoProcessor::collectData(double step)
{
   values[count++] = step;
   for (SPtr<Mic> mic : microphones)
   {
      LBMReal f[D3Q27System::ENDF+1];
      mic->distridution->getDistribution(f, val<1>(mic->nodeIndexes), val<2>(mic->nodeIndexes), val<3>(mic->nodeIndexes));
      LBMReal vx1, vx2, vx3, rho;
      calcMacros(f, rho, vx1, vx2, vx3);
      values[count++] = rho;
   }
}

void MicrophoneArrayCoProcessor::writeFile(double step)
{
   std::string fname = path+"/mic/mic_"+UbSystem::toString(comm->getProcessID())+".csv";
   std::ofstream ostr;
   ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
   if (!ostr)
   {
      ostr.clear();
      std::string path = UbSystem::getPathFromString(fname);
      if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
      if (!ostr) throw UbException(UB_EXARGS, "couldn't open file "+fname);
   }

   int line_size = static_cast<int>(microphones.size()+1);
   int total_size = static_cast<int>(values.size());

   for (int i = 0; i < total_size; i+=line_size)
   {
      int index = 0;
      for (int j = 0; j < line_size; j++)
      {
         ostr << values[i+j] << ";";
      }
      ostr << "\n";
   }

   ostr.close();
   count = 0;
   //UBLOG(logINFO,"MicrophoneArrayCoProcessor step: " << static_cast<int>(step));
}
