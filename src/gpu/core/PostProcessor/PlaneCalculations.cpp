#include "PlaneCalculations.h"

#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <cstdio>
#include <fstream>
#include <sstream>

#include <basics/StringUtilities/StringUtil.h>

void setSizeOfPlane(Parameter* para, int lev, unsigned int z)
{
   para->getParH(lev)->sizePlanePress  = 0;
   para->getParH(lev)->isSetPress      = false;
   unsigned int mm[8];
   unsigned int k = z;

   for (unsigned int j=1; j<para->getParH(lev)->gridNY + 2 * STARTOFFY - 1; j++)
   {
      for (unsigned int i=1; i<para->getParH(lev)->gridNX + 2 * STARTOFFX - 1; i++)
      {
         mm[0]= para->getParH(lev)->nx*(para->getParH(lev)->ny*k + j) + i;
         mm[1]= mm[0]                                                                       -1; //dM00
         mm[2]= mm[0]                                                -para->getParH(lev)->nx-1; //dMM0
         mm[3]= mm[0]                                                -para->getParH(lev)->nx;   //d0M0
         mm[4]= mm[0]-(para->getParH(lev)->nx*para->getParH(lev)->ny);                          //d00M
         mm[5]= mm[0]-(para->getParH(lev)->nx*para->getParH(lev)->ny)                       -1; //dM0M
         mm[6]= mm[0]-(para->getParH(lev)->nx*para->getParH(lev)->ny)-para->getParH(lev)->nx;   //d0MM
         mm[7]= mm[0]-(para->getParH(lev)->nx*para->getParH(lev)->ny)-para->getParH(lev)->nx-1; //dMMM

         if ( para->getParH(lev)->geo[mm[0]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[1]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[2]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[3]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[4]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[5]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[6]] != GEO_VOID ||
              para->getParH(lev)->geo[mm[7]] != GEO_VOID )
         {
            para->getParH(lev)->sizePlanePress  += 1;
            if (para->getParH(lev)->isSetPress == false)
            {
               para->getParH(lev)->startP = mm[0];
               para->getParH(lev)->isSetPress = true;
            }
         }
      }
   }
}


void calcPressure(Parameter* para, std::string inorout, int lev)
{
    unsigned int m   = para->getParH(lev)->startP;
   unsigned int anz = 0;
    double sumrho = 0.0, mrho = 0.0;
   double PressIn, PressOut;

   for (unsigned int i = 0; i < para->getParH(lev)->sizePlanePress; i++)
   {
      sumrho += para->getParH(lev)->rho[m];
      anz++;
      m++;
   }
    mrho = sumrho/anz;
    double deltav = para->getVelocityRatio();
    std::cout << "\n sumrho = "   << sumrho ;
    std::cout << "\n anz = "      << anz ;
    std::cout << "\n mrho = "     << mrho ;
    std::cout << "\n deltarho = " << para->getDensityRatio();
    std::cout << "\n deltav = "   << deltav ;
    if (inorout=="in")
    {
        PressIn = para->getDensityRatio() * deltav * deltav * mrho/3.0;
        std::cout << "\n Druck Ein: " << PressIn << " Pa \n";
    }
    else if (inorout=="out")
    {
        PressOut = para->getDensityRatio() * deltav * deltav * mrho/3.0;
        std::cout << "\n Druck Aus: " << PressOut << " Pa \n";
    }
}




void calcFlowRate(Parameter* para, int lev)
{
   unsigned int m, sizePlane;
   m = para->getParH(lev)->startPIN;
   sizePlane = para->getParH(lev)->sizePlanePressIN;

   unsigned int anz = 0;
   double FlowRate = 0.0;
   double sumvelo = 0.0, mvelo = 0.0;

   for (unsigned int i = 0; i < sizePlane; i++)
   {
      if (para->getParH(lev)->geo[m] == GEO_FLUID)
      {
         sumvelo += para->getParH(lev)->velocityZ[para->getParH(lev)->k[m]];
         anz++;
      } 
      m++;
   }
   mvelo = sumvelo/anz;
   //double deltav = ic.u0_ratio;
   double RealVelM = mvelo * sqrt(3.0) * 349.08 /*cs*/; //[m/s]
   double RealX = (double)para->getRealX(); //[m]
   double RealY = (double)para->getRealY(); //[m]
   std::cout << "\n sizePlane = " << sizePlane ;
   std::cout << "\n sumvelo = "   << sumvelo ;
   std::cout << "\n anz = "       << anz ;
   std::cout << "\n mvelo = "     << mvelo ;
   std::cout << "\n RealVelM = "  << RealVelM ;
   FlowRate = RealVelM * RealX * RealY * 3600./*s->h*/;
   std::cout << "\n Flow Rate: " << FlowRate << " m^3/h \n";
}




































//////////////////////////////////////////////////////////////////////////
//advection + diffusion
//////////////////////////////////////////////////////////////////////////

void calcPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager, int lev)
{
    //////////////////////////////////////////////////////////////////////////
    //copy to host
    //coarsest Grid ... with the pressure nodes
    //please test -> Copy == Alloc ??
    ////////////////////////////////////////////
    //Version Press neighbor
    unsigned int NoNin   = para->getParH(lev)->numberOfPointsCpTop;
    unsigned int NoNout1 = para->getParH(lev)->numberOfPointsCpBottom;
    unsigned int NoNout2 = para->getParH(lev)->pressureBC.numberOfBCnodes;
    ////////////////////////////////////////////
    ////Version cp top
    //unsigned int NoN = para->getParH(lev)->numberOfPointsCpTop;
    ////////////////////////////////////////////
    ////Version cp bottom
    //unsigned int NoN = para->getParH(lev)->numberOfPointsCpBottom;

    cudaMemoryManager->cudaCopyPlaneConcIn(lev, NoNin);
    cudaMemoryManager->cudaCopyPlaneConcOut1(lev, NoNout1);
    cudaMemoryManager->cudaCopyPlaneConcOut2(lev, NoNout2);
    ////////////////////////////////////////////
    //calculate concentration
    double concPlaneIn = 0.;
    double concPlaneOut1 = 0.;
    double concPlaneOut2 = 0.;
    ////////////////////////////////////////////
    double counter1 = 0.;
    for (unsigned int it = 0; it < NoNin; it++)
    {
        if (para->getParH(lev)->typeOfGridNode[it] == GEO_FLUID)
        {
            concPlaneIn   += (double) (para->getParH(lev)->ConcPlaneIn[it]);
            counter1 += 1.;
        }
    }
    concPlaneIn /= (double)(counter1);
    ////////////////////////////////////////////
    counter1 = 0.;
    for (unsigned int it = 0; it < NoNout1; it++)
    {
        if (para->getParH(lev)->typeOfGridNode[it] == GEO_FLUID)
        {
            concPlaneOut1 += (double) (para->getParH(lev)->ConcPlaneOut1[it]);
            counter1 += 1.;
        }
    }
    concPlaneOut1 /= (double)(counter1);
    ////////////////////////////////////////////
    counter1 = 0.;
    for (unsigned int it = 0; it < NoNout2; it++)
    {
        if (para->getParH(lev)->typeOfGridNode[it] == GEO_FLUID)
        {
            concPlaneOut2 += (double) (para->getParH(lev)->ConcPlaneOut2[it]);
            counter1 += 1.;
        }
    }
    concPlaneOut2 /= (double)(counter1);
    ////////////////////////////////////////////
    //concPlaneIn /= (double)(NoN);
    //concPlaneOut1 /= (double)(NoN);
    //concPlaneOut2 /= (double)(NoN);
    //////////////////////////////////////////////////////////////////////////
    //Copy to vector x,y,z
    para->getParH(lev)->PlaneConcVectorIn.push_back(concPlaneIn);
    para->getParH(lev)->PlaneConcVectorOut1.push_back(concPlaneOut1);
    para->getParH(lev)->PlaneConcVectorOut2.push_back(concPlaneOut2);
    //////////////////////////////////////////////////////////////////////////
}



void allocPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    //////////////////////////////////////////////////////////////////////////
    //set level   ---> maybe we need a loop
    int lev = para->getCoarse();
    //////////////////////////////////////////////////////////////////////////
    //allocation
    //coarsest Grid ... with the pressure nodes
    //please test -> Copy == Alloc ??
    ////////////////////////////////////////////
    //Version Press neighbor
    cudaMemoryManager->cudaAllocPlaneConcIn(lev, para->getParH(lev)->numberOfPointsCpTop);
    cudaMemoryManager->cudaAllocPlaneConcOut1(lev, para->getParH(lev)->numberOfPointsCpBottom);
    cudaMemoryManager->cudaAllocPlaneConcOut2(lev, para->getParH(lev)->pressureBC.numberOfBCnodes);
    printf("\n Number of elements plane concentration = %d + %d + %d \n", para->getParH(lev)->numberOfPointsCpTop, para->getParH(lev)->numberOfPointsCpBottom, para->getParH(lev)->pressureBC.numberOfBCnodes);
    ////////////////////////////////////////////
    ////Version cp top
    //para->cudaAllocPlaneConc(lev, para->getParH(lev)->numberOfPointsCpTop);
    //printf("\n Number of elements plane concentration = %d \n", para->getParH(lev)->numberOfPointsCpTop);
    ////////////////////////////////////////////
    ////Version cp bottom
    //para->cudaAllocPlaneConc(lev, para->getParH(lev)->numberOfPointsCpBottom);
    //printf("\n Number of elements plane concentration = %d \n", para->getParH(lev)->numberOfPointsCpBottom);
    //////////////////////////////////////////////////////////////////////////
}



void printPlaneConc(Parameter* para, CudaMemoryManager* cudaMemoryManager)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //set level   ---> maybe we need a loop
    int lev = para->getCoarse();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffnameIn = para->getFName() + UbSystem::toString(para->getMyProcessID()) + "_" + "In" + "_PlaneConc.txt";
    const char* fnameIn = ffnameIn.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostrIn;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostrIn.open(fnameIn);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    for (size_t i = 0; i < para->getParH(lev)->PlaneConcVectorIn.size(); i++)
    {
        ostrIn << para->getParH(lev)->PlaneConcVectorIn[i]  << std::endl ;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostrIn.close();
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffnameOut1 = para->getFName() + UbSystem::toString(para->getMyProcessID()) + "_" + "Out1" + "_PlaneConc.txt";
    const char* fnameOut1 = ffnameOut1.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostrOut1;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostrOut1.open(fnameOut1);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    for (size_t i = 0; i < para->getParH(lev)->PlaneConcVectorOut1.size(); i++)
    {
        ostrOut1 << para->getParH(lev)->PlaneConcVectorOut1[i]  << std::endl ;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostrOut1.close();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffnameOut2 = para->getFName() + UbSystem::toString(para->getMyProcessID()) + "_" + "Out2" + "_PlaneConc.txt";
    const char* fnameOut2 = ffnameOut2.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostrOut2;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostrOut2.open(fnameOut2);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    for (size_t i = 0; i < para->getParH(lev)->PlaneConcVectorOut2.size(); i++)
    {
        ostrOut2 << para->getParH(lev)->PlaneConcVectorOut2[i]  << std::endl ;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostrOut2.close();
    //////////////////////////////////////////////////////////////////////////
    cudaMemoryManager->cudaFreePlaneConc(lev);
    //////////////////////////////////////////////////////////////////////////
}





//////////////////////////////////////////////////////////////////////////
//Print Test round of Error
void printRE(Parameter* para, CudaMemoryManager* cudaMemoryManager, int timestep)
{
    //////////////////////////////////////////////////////////////////////////
    //set level
    int lev = 0;
    //////////////////////////////////////////////////////////////////////////
    //set filename
    std::string ffname = para->getFName()+StringUtil::toString<int>(para->getMyProcessID())+"_"+StringUtil::toString<int>(timestep)+"_RE.txt";
    const char* fname = ffname.c_str();
    //////////////////////////////////////////////////////////////////////////
    //set ofstream
    std::ofstream ostr;
    //////////////////////////////////////////////////////////////////////////
    //open file
    ostr.open(fname);
    //////////////////////////////////////////////////////////////////////////
    //fill file with data
    bool doNothing = false;
    for (unsigned int i = 0; i < para->getParH(lev)->pressureBC.numberOfBCnodes; i++)
    {
        doNothing = false;
        for (std::size_t j = 0; j < 27; j++)
        {
            if (para->getParH(lev)->kDistTestRE.f[0][j*para->getParH(lev)->pressureBC.numberOfBCnodes + i]==0)
            {
                doNothing = true;
                continue;
            }
            ostr << para->getParH(lev)->kDistTestRE.f[0][j*para->getParH(lev)->pressureBC.numberOfBCnodes + i]  << "\t";
        }
        if (doNothing==true)
        {
            continue;
        }
        ostr << std::endl;
    }
    //////////////////////////////////////////////////////////////////////////
    //close file
    ostr.close();
    //////////////////////////////////////////////////////////////////////////
    if (timestep == (int)para->getTimestepEnd())
    {
        cudaMemoryManager->cudaFreeTestRE(lev);
    }
    //////////////////////////////////////////////////////////////////////////
}
