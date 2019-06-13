#include "Calculation/DragLift.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include <stdio.h>

#include <fstream>
#include <sstream>

//#include <math.h>
//#include "LB.h"

using namespace std;



void calcDragLift(Parameter* para, CudaMemoryManager* cudaManager, int lev)
{
	//////////////////////////////////////////////////////////////////////////
	//copy to host
	//finest Grid ... with the geometry nodes
	//please test -> Copy == Alloc ??
	cudaManager->cudaCopyDragLift(lev, para->getParH(lev)->QGeom.kQ);
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
	double delta_x_F = 0.00703125;//[m] fine  16.22MI
	//double delta_x_F = 0.00625;//[m] fine  16.22MI
	//double delta_x_F = 0.0046875;//[m] fine  16.43MI
	//double delta_x_F = 0.003125;//[m] fine  16.117MI
	//////////////////////////////////////////////////////////////////////////
	double A  = 2.16693/(delta_x_F*delta_x_F);// Car 
	//////////////////////////////////////////////////////////////////////////

	//double LBtoSI = 1.0;//Sphere 
	//double A  = 110.0 * 28.0; //Ship width times height in fine nodes
	//double delta_x = 0.0045;//[m] fine
	//double delta_t = para->getVelocity() * delta_x / 15.96; 
	//double LBtoSI = 1.204 * (pow(delta_x, 4))/(pow(delta_t,2));//rho_SI * delta_x^4 / delta_t^2 = 1.204 kg/m³ * (0.0045m)^4 / (0.00000757s)^2 ... LB to kg*m/s²
	//double LBtoSI = 1000 * (pow(delta_x, 4))/(pow(delta_t,2));//rho_SI * delta_x^4 / delta_t^2 = 1000 kg/m³ * (0.1m)^4 / (0.00187s)^2 ... LB to kg*m/s²

	for (int it = 0; it < para->getParH(lev)->QGeom.kQ; it++)
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



void allocDragLift(Parameter* para, CudaMemoryManager* cudaManager)
{
	//////////////////////////////////////////////////////////////////////////
	//set level
	int lev = para->getMaxLevel();
	//////////////////////////////////////////////////////////////////////////
	//allocation
	//finest Grid ... with the geometry nodes
	//please test -> Copy == Alloc ??
	cudaManager->cudaAllocDragLift(lev, para->getParH(lev)->QGeom.kQ);
	//////////////////////////////////////////////////////////////////////////
	printf("\n Anzahl Elemente fuer Drag Lift = %d \n", para->getParH(lev)->QGeom.kQ);
}



void printDragLift(Parameter* para, CudaMemoryManager* cudaManager, int timestep)
{
	//////////////////////////////////////////////////////////////////////////
	//set level
	int lev = para->getMaxLevel();
	//////////////////////////////////////////////////////////////////////////
	//set filename
	std::string ffname = para->getFName()+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(timestep)+"_DragLift.txt";
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
	if (timestep == para->getTEnd())
	{
		cudaManager->cudaFreeDragLift(lev);
	}
	//////////////////////////////////////////////////////////////////////////
}