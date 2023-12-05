#include "ForceCalculations.h"

//////////////////////////////////////////////////////////////////////////
#include "GPU/GPU_Interface.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "GPU/CudaMemoryManager.h"

#include "PostProcessor/MacroscopicQuantities.cuh"

#include "StringUtilities/StringUtil.h"
//using namespace std;
//////////////////////////////////////////////////////////////////////////

ForceCalculations::ForceCalculations(Parameter* para)
{
    Ta = 10.0;        // number of time steps between adjusting
    vx1Targed = para->getVelocity(); // objected LB velocity

    Kpcrit = 3.0 / Ta;// 0.3;
    Tcrit = 3.0 * Ta; // 30.0;
    Tn = 1.0 * Tcrit; //0.5 * Tcrit;
    Tv = 0.24 * Tcrit; //0.12 * Tcrit;

    Kp = 0.6 * Kpcrit;
    Ki = Kp / Tn;
    Kd = Kp * Tv;

    y = 0.0;
    e = 0.0;
    esum = 0.0;
    eold = 0.0;

    isPID = true;
}

void ForceCalculations::calcPIDControllerForForce(Parameter* para, CudaMemoryManager* cudaMemoryManager)
 {
     //////////////////////////////////////////////////////////////////////////
     double tempVeloX = 0.0;
     double veloAverageX = 0.0;
     double levelVeloAverageX = 0.0;
     int counter = 0;
     //////////////////////////////////////////////////////////////////////////
     for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
     {
         //////////////////////////////////////////////////////////////////////
         //measure the velocity
         unsigned long long numberOfElements = para->getParH(lev)->numberOfNodes;
         if (numberOfElements > 0)
         {
             CalcMacCompSP27(para->getParD(lev)->velocityX,
                             para->getParD(lev)->velocityY,
                             para->getParD(lev)->velocityZ,
                             para->getParD(lev)->rho,
                             para->getParD(lev)->pressure,
                             para->getParD(lev)->typeOfGridNode,
                             para->getParD(lev)->neighborX,
                             para->getParD(lev)->neighborY,
                             para->getParD(lev)->neighborZ,
                             para->getParD(lev)->numberOfNodes,
                             para->getParD(lev)->numberofthreads,
                             para->getParD(lev)->distributions.f[0],
                             para->getParD(lev)->isEvenTimestep);
             getLastCudaError("CalcMacSP27 execution failed");
             //////////////////////////////////////////////////////////////////
             cudaMemoryManager->cudaCopyPrint(lev);
             //////////////////////////////////////////////////////////////////
             for (size_t pos = 0; pos < numberOfElements; pos++)
             {
                 tempVeloX += (double)para->getParH(lev)->velocityX[pos];
             }
             tempVeloX /= (double)numberOfElements;
             //////////////////////////////////////////////////////////////////
             levelVeloAverageX += tempVeloX;
             //////////////////////////////////////////////////////////////////
             counter++;
             //////////////////////////////////////////////////////////////////
         }
     }
     //////////////////////////////////////////////////////////////////////////
     veloAverageX = levelVeloAverageX / (double)counter;
     //////////////////////////////////////////////////////////////////////////
     if (isPID)
     {
         //PID-Controller
         e = vx1Targed - veloAverageX;
         esum = esum + e;
         y = Kp * e + Ki * Ta * esum + Kd * (e - eold) / Ta;
         eold = e;

         y = y / 2.0;
     }
     //////////////////////////////////////////////////////////////////////////
     para->getForcesDouble()[0] = (para->getForcesDouble()[0] + y);
     para->getForcesDouble()[1] = (para->getForcesDouble()[1] + y) * 0.0;
     para->getForcesDouble()[2] = (para->getForcesDouble()[2] + y) * 0.0;
     //////////////////////////////////////////////////////////////////////////
     para->getForcesHost()[0] = (real)(para->getForcesHost()[0] + y);
     para->getForcesHost()[1] = (real)(para->getForcesHost()[1] + y) * (real)0.0;
     para->getForcesHost()[2] = (real)(para->getForcesHost()[2] + y) * (real)0.0;
     //////////////////////////////////////////////////////////////////////////
     cudaMemoryManager->cudaCopyForcingToDevice();
     //////////////////////////////////////////////////////////////////////////
 }


void ForceCalculations::printForcing(Parameter* para)
{
    //////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffname = para->getFName() + StringUtil::toString<int>(para->getMyProcessID()) + "_forcing.txt";
    const char* fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostr.open(fname, std::fstream::app);
    ostr << para->getForcesHost()[0] << " " << para->getForcesHost()[1] << " " << para->getForcesHost()[2];
    ostr << std::endl;
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
}
