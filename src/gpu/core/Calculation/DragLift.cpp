//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================

#include "DragLift.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cstdio>
#include <fstream>
#include <sstream>

#include "GPU/CudaMemoryManager.h"
#include "Parameter/Parameter.h"
#include "StringUtilities/StringUtil.h"

using namespace std;

void calcDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev)
{
    //////////////////////////////////////////////////////////////////////////
    //copy to host
    //finest Grid ... with the geometry nodes
    //please test -> Copy == Alloc ??
    cudaMemoryManager->cudaCopyDragLift(lev, para->getParH(lev)->geometryBC.numberOfBCnodes);
    //////////////////////////////////////////////////////////////////////////
    //calc drag
    double dragX = 0., dragY = 0., dragZ = 0.;
    double CDX   = 0., CDY   = 0., CDZ   = 0.;
    //double Pi = 3.14159265358979323846;
    //double A  = Pi * 1000.0 * 1000.0;// Sphere 
    //////////////////////////////////////////////////////////////////////////
    //Cars
    //double delta_x_F = 0.00625;//[m] fine  11.25MI
    //double delta_x_F = 0.0045;//[m] fine  12.54MI
    //double delta_x_F = 0.00359375;//[m] fine  13.54MI
    ///////////////////////////////
    //double delta_x_F = 0.00703125;//[m] fine  16.22MI
    ////double delta_x_F = 0.00625;//[m] fine  16.22MI
    ////double delta_x_F = 0.0046875;//[m] fine  16.43MI
    ////double delta_x_F = 0.003125;//[m] fine  16.117MI
    ///////////////////////////////
    //DLC
    double delta_x_F = 0.00625;//[m]

    //////////////////////////////////////////////////////////////////////////
    //double A  = 2.16693/(delta_x_F*delta_x_F);// Car 
    double A  = 2.19/(delta_x_F*delta_x_F);// DLC 
    //////////////////////////////////////////////////////////////////////////

    //double LBtoSI = 1.0;//Sphere 
    //double A  = 110.0 * 28.0; //Ship width times height in fine nodes
    //double delta_x = 0.0045;//[m] fine
    //double delta_t = para->getVelocity() * delta_x / 15.96; 
    //double LBtoSI = 1.204 * (pow(delta_x, 4))/(pow(delta_t,2));//rho_SI * delta_x^4 / delta_t^2 = 1.204 kg/m^3 * (0.0045m)^4 / (0.00000757s)^2 ... LB to kg*m/s^2
    //double LBtoSI = 1000 * (pow(delta_x, 4))/(pow(delta_t,2));//rho_SI * delta_x^4 / delta_t^2 = 1000 kg/m^3 * (0.1m)^4 / (0.00187s)^2 ... LB to kg*m/s^2

    for (unsigned int it = 0; it < para->getParH(lev)->geometryBC.numberOfBCnodes; it++)
    {
        dragX += (double) (para->getParH(lev)->DragPreX[it] - para->getParH(lev)->DragPostX[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
        dragY += (double) (para->getParH(lev)->DragPreY[it] - para->getParH(lev)->DragPostY[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
        dragZ += (double) (para->getParH(lev)->DragPreZ[it] - para->getParH(lev)->DragPostZ[it]); //Kraft da Impuls pro Zeitschritt merke: andere nennen es FD
    }
    //////////////////////////////////////////////////////////////////////////
    //calc CD
    CDX = 2.0 * dragX / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    CDY = 2.0 * dragY / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    CDZ = 2.0 * dragZ / (1.0 /*rho_0*/ * para->getVelocity() * para->getVelocity() * A);
    //////////////////////////////////////////////////////////////////////////
    //transform CD to SI
    //CDX *= LBtoSI;
    //CDY *= LBtoSI;
    //CDZ *= LBtoSI;
    //////////////////////////////////////////////////////////////////////////
    //Copy to vector x,y,z
    para->getParH(lev)->DragXvector.push_back(CDX);
    para->getParH(lev)->DragYvector.push_back(CDY);
    para->getParH(lev)->DragZvector.push_back(CDZ);
    //////////////////////////////////////////////////////////////////////////
}



void allocDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    //////////////////////////////////////////////////////////////////////////
    //set level
    int lev = para->getMaxLevel();
    //////////////////////////////////////////////////////////////////////////
    //allocation
    //finest Grid ... with the geometry nodes
    //please test -> Copy == Alloc ??
    cudaMemoryManager->cudaAllocDragLift(lev, para->getParH(lev)->geometryBC.numberOfBCnodes);
    //////////////////////////////////////////////////////////////////////////
    printf("\n Anzahl Elemente fuer Drag Lift = %d \n", para->getParH(lev)->geometryBC.numberOfBCnodes);
}



void printDragLift(Parameter* para, CudaMemoryManager* cudaMemoryManager, int timestep)
{
    //////////////////////////////////////////////////////////////////////////
    //set level
    int lev = para->getMaxLevel();
    //////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffname = para->getFName()+StringUtil::toString<int>(para->getMyProcessID())+"_"+StringUtil::toString<int>(timestep)+"_DragLift.txt";
    const char* fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostr.open(fname);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    for (size_t i = 0; i < para->getParH(lev)->DragXvector.size(); i++)
    {
        ostr << para->getParH(lev)->DragXvector[i] << "\t" << para->getParH(lev)->DragYvector[i] << "\t" << para->getParH(lev)->DragZvector[i] << endl ;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
    if (timestep == (int)para->getTimestepEnd())
    {
        cudaMemoryManager->cudaFreeDragLift(lev);
    }
    //////////////////////////////////////////////////////////////////////////
}