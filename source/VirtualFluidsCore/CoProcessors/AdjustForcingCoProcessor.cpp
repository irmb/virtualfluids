/*
* D3Q27AdjustForcingCoProcessor.cpp
*
*
*  Author: Konstantin Kutscher
*/

#include "AdjustForcingCoProcessor.h"

#include <SetForcingBlockVisitor.h>

#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

using namespace std;

AdjustForcingCoProcessor::AdjustForcingCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
   const std::string& path,
   IntegrateValuesHelperPtr integrateValues,
   double vTarged,
   CommunicatorPtr comm)

   : CoProcessor(grid, s),
   path(path),
   integrateValues(integrateValues),
   comm(comm),
   vx1Targed(vTarged),
   forcing(forcing)
{
   //cnodes = integrateValues->getCNodes();
   root = comm->isRoot();

   Ta = scheduler->getMaxStep();

   Kpcrit = 3.0 / Ta;// 0.3;
   Tcrit = 3.0 * Ta; // 30.0;
   Tn = 0.5 * Tcrit;
   Tv = 0.12 * Tcrit;

   Kp = 0.6 * Kpcrit;
   Ki = Kp / Tn;
   Kd = Kp * Tv;

   y = 0;
   e = 0;
   esum = 0;
   eold = 0;
   forcing = 0;

   if (root)
   {
      std::string fname = path + "/forcing/forcing.csv";
      std::ofstream ostr;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if (!ostr)
      {
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if (path.size() > 0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file " + fname);
      }
      ostr << "step;volume;vx1average;forcing\n";
      ostr.close();

      //////////////////////////////////////////////////////////////////////////////////////////////////
      //temporere Lösung
      std::string fNameCfg = path + "/forcing/forcing.cfg";
      std::ifstream istr2;
      istr2.open(fNameCfg.c_str(), std::ios_base::in);
      if (istr2)
      {
         istr2 >> forcing;
         //istr2 >> esum;
         //istr2 >> eold;
      }
      istr2.close();
   }
   ////////////////////////////////////////////////////////////////////////////////////////////////////////
}
//////////////////////////////////////////////////////////////////////////
AdjustForcingCoProcessor::~AdjustForcingCoProcessor()
{
}
//////////////////////////////////////////////////////////////////////////
void AdjustForcingCoProcessor::process(double step)
{
   if (scheduler->isDue(step))
      collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void AdjustForcingCoProcessor::collectData(double step)
{
   //////////////////////////////////////////////////////////////////////////////////////////////////
   //temporere Lösung
   if (root)
   {
      std::string fNameCfg = path + "/forcing/forcing.cfg";
      std::ofstream ostr2;
      ostr2.open(fNameCfg.c_str(), std::ios_base::out);
      if (!ostr2)
      {
         ostr2.clear();
         string path = UbSystem::getPathFromString(fNameCfg);
         if (path.size() > 0) { UbSystem::makeDirectory(path); ostr2.open(fNameCfg.c_str(), std::ios_base::out); }
         if (!ostr2) throw UbException(UB_EXARGS, "couldn't open file " + fNameCfg);
      }
      ostr2 << forcing << " " << esum << " " << eold;
      ostr2.close();
   }
   ////////////////////////////////////////////////////////////////////////////////////////////////////////

   integrateValues->calculateMQ();

   if (root)
   {
      cellsVolume = integrateValues->getCellsVolume();
      double vx1 = integrateValues->getVx1();
      vx1Average = (vx1 / cellsVolume);

      //////////////////////////////////////////////////////////////////////////
      //PID-Controller (PID-Regler)
      e = vx1Targed - vx1Average;
      esum = esum + e;
      y = Kp * e + Ki * Ta * esum + Kd * (e - eold) / Ta;
      eold = e;

      forcing = forcing + y;
      //////////////////////////////////////////////////////////////////////////
   }
   //////////////////////////////////////////////////////////////////////////
   comm->broadcast(forcing);

   mu::Parser fctForcingX1, fctForcingX2, fctForcingX3;
   fctForcingX1.SetExpr("Fx1");
   fctForcingX1.DefineConst("Fx1", forcing);
   fctForcingX2.SetExpr("0.0");
   fctForcingX3.SetExpr("0.0");
   SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
   grid->accept(forcingVisitor);

   //BOOST_FOREACH(CalcNodes cn, cnodes)
   //{
   //   LBMKernel3DPtr kernel = cn.block->getKernel();
   //   if (kernel)
   //   {
   //      kernel->setForcingX1(fctForcingX1);
   //      kernel->setWithForcing(true);
   //   }
   //      
   //}

   if (root)
   {
      //UBLOG(logINFO, "D3Q27AdjustForcingCoProcessor step: " << static_cast<int>(step));
      //UBLOG(logINFO, "new forcing is: " << forcing);
      std::string fname = path + "/forcing/forcing.csv";
      //std::string fname = path + "/forcing/forcing_"+UbSystem::toString(comm->getProcessID())+".csv";
      std::ofstream ostr;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if (!ostr)
      {
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if (path.size() > 0) { UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app); }
         if (!ostr) throw UbException(UB_EXARGS, "couldn't open file " + fname);
      }
      int istep = static_cast<int>(step);

      ostr << istep << ";" << cellsVolume << ";" << vx1Average << "; " << forcing << "\n";
      ostr.close();

   }
}
