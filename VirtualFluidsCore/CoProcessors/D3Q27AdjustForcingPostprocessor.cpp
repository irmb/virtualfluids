/*
* D3Q27AdjustForcingPostprocessor.cpp
*
*  
*  Author: Sonja Uphoff
*/

#include "D3Q27AdjustForcingPostprocessor.h"

#include <SetForcingBlockVisitor.h>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

using namespace std;

D3Q27AdjustForcingPostprocessor::D3Q27AdjustForcingPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                                                                 const std::string& path,
                                                                 D3Q27IntegrateValuesHelperPtr integrateValues, 
                                                                 LBMReal vTarged,
                                                                 LBMReal forcing,
                                                                 CommunicatorPtr comm)

                                                                 : Postprocessor(grid, s),
                                                                 path(path),
                                                                 integrateValues(integrateValues),
                                                                 comm(comm),
                                                                 vTarged(vTarged),
                                                                 vPreviousStep(0.0),
                                                                 forcing(forcing)
{
   cnodes = integrateValues->getCNodes();
   if (comm->getProcessID() == comm->getRoot())
   {
      std::string fname = path+"/forcing/forcing.csv";
      std::ofstream ostr;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if (!ostr)
      {
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file "+fname);
      }
      ostr << "step;volume;vx1average;factor;forcing\n";
      ostr.close();
   }
}
//////////////////////////////////////////////////////////////////////////
D3Q27AdjustForcingPostprocessor::~D3Q27AdjustForcingPostprocessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void D3Q27AdjustForcingPostprocessor::update(double step)
{
   if(scheduler->isDue(step) )
      collectPostprocessData(step);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27AdjustForcingPostprocessor::collectPostprocessData(double step)
{
   integrateValues->calculateMQ();

   UBLOG(logDEBUG3, "D3Q27AdjustForcingPostprocessor::update:" << step);
   int gridRank = grid->getRank();
   int minInitLevel = this->grid->getCoarsestInitializedLevel();
   int maxInitLevel = this->grid->getFinestInitializedLevel();

   double cellsVolume = integrateValues->getCellsVolume();

   double vx1 = integrateValues->getVx1();
   double vx1average = (vx1/cellsVolume);


   double newForcing = forcing;

   if (!((vPreviousStep>(vx1average) && (vx1average)>vTarged)  || (vPreviousStep<(vx1average) && ((vx1average)<vTarged))))
   {
      double C = 1.0; //0.7 //P; //free parameter [0.1, 1]
      newForcing = forcing*((1-((vx1average-vTarged)/vTarged)*C));
      newForcing=UbMath::max(newForcing,0.0);
      newForcing=UbMath::min(newForcing,5e-3);
   }

   vPreviousStep=vx1average;

   forcing = newForcing;

   mu::Parser fctForcingX1, fctForcingX2, fctForcingX3;
   fctForcingX1.SetExpr("Fx1");
   fctForcingX1.DefineConst("Fx1", newForcing);
   fctForcingX2.SetExpr("0.0");
   fctForcingX3.SetExpr("0.0");
   //SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
   //grid->accept(forcingVisitor);

   BOOST_FOREACH(CalcNodes cn, cnodes)
   {
      LBMKernel3DPtr kernel = cn.block->getKernel();
      if (kernel)
      {
         kernel->setForcingX1(fctForcingX1);
         kernel->setWithForcing(true);
      }
         
   }

   if (comm->getProcessID() == comm->getRoot())
   {
      //UBLOG(logINFO, "D3Q27AdjustForcingPostprocessor step: " << static_cast<int>(step));
      //UBLOG(logINFO, "new forcing is: " << forcing);
      std::string fname = path+"/forcing/forcing.csv";
      std::ofstream ostr;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if (!ostr)
      {
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if (path.size()>0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file "+fname);
      }
      int istep = static_cast<int>(step);
      ostr << istep << ";" << cellsVolume << ";" << vx1average << "; " << 1-((vx1average-vTarged)/vTarged) << "; " << forcing << "\n";
      ostr.close();
   }
}
