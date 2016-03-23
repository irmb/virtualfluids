/*
 * D3Q27RhoPostprocessor.cpp
 *
 *  Created on: 28.12.2010
 *      Author: kucher
 */

#include "D3Q27PressureDifferencePostprocessor.h"


#include <iostream>
#include <fstream>

using namespace std;

D3Q27PressureDifferencePostprocessor::D3Q27PressureDifferencePostprocessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path,
                                                                 D3Q27IntegrateValuesHelperPtr h1, D3Q27IntegrateValuesHelperPtr h2, 
                                                                 LBMReal rhoReal, LBMReal uReal, LBMReal uLB,
                                                                 CommunicatorPtr comm)

                                                : Postprocessor(grid, s)
                                                , path(path)
																, h1(h1)
																, h2(h2)
                                                ,comm(comm)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      std::ofstream ostr;
      string fname = path;
      ostr.open(fname.c_str(), std::ios_base::out);
      if(!ostr)
      { 
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }
      ostr << "step" << "\t" << "nodes1" << "\t" << "nodes2" << "\t";
      ostr << "sRho1" << "\t" << "p1_1"  << "\t" << "sRho2" << "\t" << "p1_2" << "\t" << "deltaP1"<< "\t";
      ostr << "sPress1" << "\t" << "p2_1" << "\t" << "sPress2" << "\t" << "p2_2" << "\t" << "deltaP2";
      ostr << endl;
      ostr.close();

      factor1 = (1.0/3.0)*rhoReal*(uReal/uLB)*(uReal/uLB);
      factor2 = rhoReal*(uReal/uLB)*(uReal/uLB);
   }
}
//////////////////////////////////////////////////////////////////////////
D3Q27PressureDifferencePostprocessor::~D3Q27PressureDifferencePostprocessor() 
{
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PressureDifferencePostprocessor::update(double step)
{
   if(scheduler->isDue(step) )
      collectPostprocessData(step);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27PressureDifferencePostprocessor::collectPostprocessData(double step)
{
   h1->calculateMQ();
   h2->calculateMQ();
   
   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = static_cast<int>(step);
      std::ofstream ostr;
      double nn1 = h1->getNumberOfFluidsNodes();
      double nn2 = h2->getNumberOfFluidsNodes();
      double rho1 = h1->getRho();
      double rho2 = h2->getRho();
      double p1_1 = (rho1/nn1) * factor1;
      double p1_2 = (rho2/nn2) * factor1;
      double dp1 = p1_1 - p1_2;

      //double press1 = h1->getPress();
      //double press2 = h2->getPress();
      //double p2_1 = (press1/nn1) * factor2;
      //double p2_2 = (press2/nn2) * factor2;
      //double dp2 = p2_1 - p2_2;

      string fname = path;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }

      ostr << istep << "\t" << nn1 << "\t"  << nn2 << "\t"; 
      ostr << rho1 << "\t" << p1_1 << "\t" << rho2 << "\t" << p1_2 << "\t" << dp1 << "\t";
      //ostr << press1 << "\t" << p2_1 << "\t" << press2 << "\t" << p2_2 << "\t" << dp2;
      ostr << endl;
      ostr.close();
   }
}
