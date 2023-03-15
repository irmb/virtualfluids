#include "CalculateTorqueCoProcessor.h"
#include "BCProcessor.h"

#include <mpi/Communicator.h>
#include "D3Q27Interactor.h"
#include "UbScheduler.h"
#include "Grid3D.h"
#include "BoundaryConditions.h"
#include "DataSet3D.h"
#include "Block3D.h"
#include "LBMKernel.h"
#include "BCArray3D.h"
#include "EsoTwist3D.h"
#include "DistributionArray3D.h"
#include "Rheology.h"

CalculateTorqueCoProcessor::CalculateTorqueCoProcessor( SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path_, std::shared_ptr<vf::mpi::Communicator> comm) : CoProcessor(grid, s), path(path_), comm(comm), torqueX1global(0), torqueX2global(0), torqueX3global(0)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      std::ofstream ostr;
       std::string fname = path_;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         const std::string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }

      ostr << "step;";
      ostr << "Tx;";
      ostr << "Ty;";
      ostr << "Tz" << std::endl;
      ostr.close();
   }
}
//////////////////////////////////////////////////////////////////////////
CalculateTorqueCoProcessor::~CalculateTorqueCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void CalculateTorqueCoProcessor::process( real step )
{
   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "D3Q27ForcesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void CalculateTorqueCoProcessor::collectData( real step )
{
   calculateForces();

   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = static_cast<int>(step);
      std::ofstream ostr;
      std::string fname = path;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         std::string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }

      ostr << istep << ";";
      ostr << torqueX1global << ";";
      ostr << torqueX2global << ";";
      ostr << torqueX3global << ";";
      ostr << Fx << ";";
      ostr << Fy << ";";
      ostr << Fz;
      ostr << std::endl;
      ostr.close();
   }
}
//////////////////////////////////////////////////////////////////////////
void CalculateTorqueCoProcessor::calculateForces()
{
   torqueX1global = 0.0;
   torqueX2global = 0.0;
   torqueX3global = 0.0;

   for(SPtr<D3Q27Interactor> interactor : interactors)
   {
      real x1Centre = interactor->getGbObject3D()->getX1Centroid();
      real x2Centre = interactor->getGbObject3D()->getX2Centroid();
      real x3Centre = interactor->getGbObject3D()->getX3Centroid();

      for(BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap())
      {
         real torqueX1 = 0.0;
         real torqueX2 = 0.0;
         real torqueX3 = 0.0;

         SPtr<Block3D> block = t.first;
         std::set< std::vector<int> >& transNodeIndicesSet = t.second;

         real deltaX = grid->getDeltaX(block);

         SPtr<ILBMKernel> kernel = block->getKernel();

         SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();          
         SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions(); 

         int ghostLayerWidth = kernel->getGhostLayerWidth();
         int minX1 = ghostLayerWidth;
         int maxX1 = (int)bcArray->getNX1() - 1 - ghostLayerWidth;
         int minX2 = ghostLayerWidth;
         int maxX2 = (int)bcArray->getNX2() - 1 - ghostLayerWidth;
         int minX3 = ghostLayerWidth;
         int maxX3 = (int)bcArray->getNX3() - 1 - ghostLayerWidth;

         for(std::vector<int> node : transNodeIndicesSet)
         {
            int x1 = node[0];
            int x2 = node[1];
            int x3 = node[2];

            //without ghost nodes
            if (x1 < minX1 || x1 > maxX1 || x2 < minX2 || x2 > maxX2 ||x3 < minX3 || x3 > maxX3 ) continue;

            if(bcArray->isFluid(x1,x2,x3)) //es kann sein, dass der node von einem anderen interactor z.B. als solid gemarkt wurde!!!
            {
               SPtr<BoundaryConditions> bc = bcArray->getBC(x1,x2,x3);
               
               Vector3D worldCoordinates = grid->getNodeCoordinates(block, x1, x2, x3);
               real rx = (worldCoordinates[0] - x1Centre) / deltaX;
               real ry = (worldCoordinates[1] - x2Centre) / deltaX;
               real rz = (worldCoordinates[2] - x3Centre) / deltaX;

               // real nx = rx / sqrt(rx * rx + ry * ry + rz * rz);
               // real ny = ry / sqrt(rx * rx + ry * ry + rz * rz);
               // real nz = rz / sqrt(rx * rx + ry * ry + rz * rz);

               UbTupleDouble3 forceVec = getForces(x1, x2, x3, distributions, bc);
               //UbTupleDouble3 forceVec = getForcesFromMoments(x1, x2, x3, kernel, distributions, bc, nx, ny, nz);
               //UbTupleDouble3 forceVec = getForcesFromStressTensor(x1, x2, x3, kernel, distributions, bc, nx, ny, nz);
               /*real*/ Fx                   = val<1>(forceVec);
               /*real*/ Fy                   = val<2>(forceVec);
               /*real*/ Fz                   = val<3>(forceVec);
              


               torqueX1 += ry * Fz - rz * Fy;
               torqueX2 += rz * Fx - rx * Fz;
               torqueX3 += rx * Fy - ry * Fx;
            }
         }

         torqueX1global += torqueX1;
         torqueX2global += torqueX2;
         torqueX3global += torqueX3;
      }
   }
   std::vector<real> values;
   std::vector<real> rvalues;
   values.push_back(torqueX1global);
   values.push_back(torqueX2global);
   values.push_back(torqueX3global);

   rvalues = comm->gather(values);
   if (comm->getProcessID() == comm->getRoot())
   {
      torqueX1global = 0.0;
      torqueX2global = 0.0;
      torqueX3global = 0.0;
      
      for (int i = 0; i < (int)rvalues.size(); i+=3)
      {
         torqueX1global += rvalues[i];
         torqueX2global += rvalues[i+1];
         torqueX3global += rvalues[i+2];
      }
   }
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 CalculateTorqueCoProcessor::getForces(int x1, int x2, int x3,  SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc)
{
   UbTupleDouble3 force(0.0,0.0,0.0);

   if(bc)
   {
      //references to tuple "force"
      real& forceX1 = val<1>(force);
      real& forceX2 = val<2>(force);
      real& forceX3 = val<3>(force);
      real f,  fnbr;

      dynamicPointerCast<EsoTwist3D>(distributions)->swap();

      for(int fdir=D3Q27System::FSTARTDIR; fdir<=D3Q27System::FENDDIR; fdir++)
      {
         if(bc->hasNoSlipBoundaryFlag(fdir) || bc->hasVelocityBoundaryFlag(fdir))
         {
            const int invDir = D3Q27System::INVDIR[fdir];
            f = dynamicPointerCast<EsoTwist3D>(distributions)->getDistributionInvForDirection(x1, x2, x3, invDir);
            fnbr = dynamicPointerCast<EsoTwist3D>(distributions)->getDistributionInvForDirection(x1+D3Q27System::DX1[invDir], x2+D3Q27System::DX2[invDir], x3+D3Q27System::DX3[invDir], fdir);

            forceX1 += (f + fnbr) * D3Q27System::DX1[invDir];
            forceX2 += (f + fnbr) * D3Q27System::DX2[invDir];
            forceX3 += (f + fnbr) * D3Q27System::DX3[invDir];
         }
      }

      dynamicPointerCast<EsoTwist3D>(distributions)->swap();
   }

   return force;
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 CalculateTorqueCoProcessor::getForcesFromMoments(int x1, int x2, int x3, SPtr<ILBMKernel> kernel, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc, real nx, real ny, real nz)
{
   using namespace vf::lbm::constant;
   UbTupleDouble3 force(0.0, 0.0, 0.0);


   if (bc) {
      real f[D3Q27System::ENDF + 1];
      distributions->getDistribution(f, x1, x2, x3);
      real collFactor = kernel->getCollisionFactor();
      real shearRate = D3Q27System::getShearRate(f, collFactor);
      real rho = D3Q27System::getDensity(f);
      real omega = Rheology::getBinghamCollFactor(collFactor, shearRate, rho);
      std::array<real, 6> moments = D3Q27System::getSecondMoments(f, omega);

      // references to tuple "force"
      real &forceX1 = val<1>(force);
      real &forceX2 = val<2>(force);
      real &forceX3 = val<3>(force);

      real mxx = (moments[0] + moments[1] + moments[2])*c1o3;
      real myy = (-c2o1 * moments[1] + moments[2] + moments[0]) * c1o3; 
      real mzz = (-c2o1 * moments[2] + moments[1] + moments[0]) * c1o3;
      real mxy = moments[3];
      real mxz = moments[4];
      real myz = moments[5];
      
      forceX1 = mxx *nx + mxy*ny + mxz*nz;
      forceX2 = mxy *nx + myy*ny + myz*nz;
      forceX3 = mxz *nx + myz*ny + mzz*nz;
   }

   return force;
}
//////////////////////////////////////////////////////////////////////////
UbTupleDouble3 CalculateTorqueCoProcessor::getForcesFromStressTensor(int x1, int x2, int x3, SPtr<ILBMKernel> kernel, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc, real nx, real ny, real nz)
{
   using namespace vf::lbm::constant;
   UbTupleDouble3 force(0.0, 0.0, 0.0);

   if (bc) {
      real f[D3Q27System::ENDF + 1];
      distributions->getDistribution(f, x1, x2, x3);
      real collFactor = kernel->getCollisionFactor();
      real shearRate = D3Q27System::getShearRate(f, collFactor);
      real rho = D3Q27System::getDensity(f);
      real omega = Rheology::getBinghamCollFactor(collFactor, shearRate, rho);
      std::array<real, 6> stress = D3Q27System::getStressTensor(f, omega);

      // references to tuple "force"
      real &forceX1 = val<1>(force);
      real &forceX2 = val<2>(force);
      real &forceX3 = val<3>(force);

      real &tauXX = stress[0];
      real &tauYY = stress[1];
      real &tauZZ = stress[2];
      real &tauXY = stress[3];
      real &tauXZ = stress[4];
      real &tauYZ = stress[5];

      //https: // journals.aps.org/pre/pdf/10.1103/PhysRevE.88.013303

      forceX1 = tauXX * nx + tauXY * ny + tauXZ * nz;
      forceX2 = tauXY * nx + tauYY * ny + tauYZ * nz;
      forceX3 = tauXZ * nx + tauYZ * ny + tauZZ * nz;
   }

   return force;
}
//////////////////////////////////////////////////////////////////////////
void CalculateTorqueCoProcessor::addInteractor( SPtr<D3Q27Interactor> interactor )
{
   interactors.push_back(interactor);
}



