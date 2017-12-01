#include "Calculation/ForceCalculations.h"

//////////////////////////////////////////////////////////////////////////
#include "GPU/GPU_Interface.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
using namespace std;
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

ForceCalculations::~ForceCalculations()
{
}



void ForceCalculations::calcPIDControllerForForce(Parameter* para)
 {
	 //////////////////////////////////////////////////////////////////////////
	 double tempVeloX = 0.0, tempVeloY = 0.0, tempVeloZ = 0.0;
	 double veloAverageX = 0.0, veloAverageY = 0.0, veloAverageZ = 0.0;
	 double levelVeloAverageX = 0.0, levelVeloAverageY = 0.0, levelVeloAverageZ = 0.0;
	 int counter = 0;
	 //////////////////////////////////////////////////////////////////////////
	 for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	 {
		 //////////////////////////////////////////////////////////////////////
		 //measure the velocity
		 int numberOfElements = para->getParH(lev)->size_Mat_SP;
		 if (numberOfElements > 0)
		 {
			 CalcMacCompSP27(para->getParD(lev)->vx_SP,
							 para->getParD(lev)->vy_SP,
							 para->getParD(lev)->vz_SP,
							 para->getParD(lev)->rho_SP,
							 para->getParD(lev)->press_SP,
							 para->getParD(lev)->geoSP,
							 para->getParD(lev)->neighborX_SP,
							 para->getParD(lev)->neighborY_SP,
							 para->getParD(lev)->neighborZ_SP,
							 para->getParD(lev)->size_Mat_SP,
							 para->getParD(lev)->numberofthreads,
							 para->getParD(lev)->d0SP.f[0],
							 para->getParD(lev)->evenOrOdd);
			 getLastCudaError("CalcMacSP27 execution failed");
			 //////////////////////////////////////////////////////////////////
			 para->cudaCopyPrint(lev);
//			 para->cudaCopyForceVelo(i,numberOfElements);
			 //////////////////////////////////////////////////////////////////
			 for (int j = 0; j < numberOfElements; j++)
			 {
				 tempVeloX += (double)para->getParH(lev)->vx_SP[j];
				 tempVeloY += (double)para->getParH(lev)->vy_SP[j];
				 tempVeloZ += (double)para->getParH(lev)->vz_SP[j];
			 }
			 tempVeloX /= (double)numberOfElements;
			 tempVeloY /= (double)numberOfElements;
			 tempVeloZ /= (double)numberOfElements;
			 //////////////////////////////////////////////////////////////////
			 levelVeloAverageX += tempVeloX;
			 levelVeloAverageY += tempVeloY;
			 levelVeloAverageZ += tempVeloZ;
			 //////////////////////////////////////////////////////////////////
			 counter++;
			 //////////////////////////////////////////////////////////////////
		 }
	 }
	 //////////////////////////////////////////////////////////////////////////
	 veloAverageX = levelVeloAverageX / (double)counter;
	 veloAverageY = levelVeloAverageY / (double)counter;
	 veloAverageZ = levelVeloAverageZ / (double)counter;
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
	 para->getForcesHost()[0] = (doubflo)(para->getForcesHost()[0] + y);
	 para->getForcesHost()[1] = (doubflo)(para->getForcesHost()[1] + y) * (doubflo)0.0;
	 para->getForcesHost()[2] = (doubflo)(para->getForcesHost()[2] + y) * (doubflo)0.0;
	 //////////////////////////////////////////////////////////////////////////
	 para->cudaCopyForcingToDevice();
	 //////////////////////////////////////////////////////////////////////////
 }
