#include "Calculation/PorousMedia.h"

//////////////////////////////////////////////////////////////////////////
//#include "GPU/GPU_Interface.h"
//#include <cuda_runtime.h>
//#include <helper_cuda.h>
//#include <stdio.h>
//#include <fstream>
//#include <sstream>
//using namespace std;
//////////////////////////////////////////////////////////////////////////

PorousMedia::PorousMedia()
{
	porosity = 0.0;
	geoID = 0;
	setLBMvaluesToZero();
	setSIvaluesToZero();
	setCoordinatesToZero();
}

PorousMedia::PorousMedia(double porosity, unsigned int geoID, double darcySI, double forchheimerSI, double dxLBM, double dtLBM)
{
	this->porosity = porosity;
	this->geoID = geoID;
	this->darcySI = darcySI;
	this->forchheimerSI = forchheimerSI;
	this->dxLBM = dxLBM;
	this->dtLBM = dtLBM;
	this->forchheimerLBM = this->forchheimerSI * this->dtLBM;
	this->darcyLBM = this->darcySI * this->dxLBM * this->dxLBM;
	setCoordinatesToZero();
}

PorousMedia::~PorousMedia()
{
}


//setter
void PorousMedia::setPorosity(double porosity) {			this->porosity = porosity; }
void PorousMedia::setDarcySI(double darcySI) {				this->darcySI = darcySI; }
void PorousMedia::setForchheimerSI(double forchheimerSI) {	this->forchheimerSI = forchheimerSI; }

void PorousMedia::setStartCoordinates(double startCoordX, double startCoordY, double startCoordZ) 
{
	this->startCoordX = startCoordX;
	this->startCoordY = startCoordY;
	this->startCoordZ = startCoordZ;
}

void PorousMedia::setEndCoordinates(double endCoordX, double endCoordY, double endCoordZ)
{
	this->endCoordX = endCoordX;
	this->endCoordY = endCoordY;
	this->endCoordZ = endCoordZ;
}

void PorousMedia::setHostNodeIDsPM(unsigned int* hostNodeIDsPM)
{
	this->hostNodeIDsPM = hostNodeIDsPM;
}

void PorousMedia::setDeviceNodeIDsPM(unsigned int* deviceNodeIDsPM)
{
	this->deviceNodeIDsPM = deviceNodeIDsPM;
}

//void PorousMedia::definePMarea(Parameter* para, unsigned int level)
//{
//	unsigned int counter = 0;
//	for (unsigned int i = 0; i < para->getParH(level)->size_Mat_SP; i++)
//	{
//		if (((para->getParH(level)->coordX_SP[i] >= this->startCoordX) && (para->getParH(level)->coordX_SP[i] <= this->endCoordX)) &&
//			((para->getParH(level)->coordY_SP[i] >= this->startCoordY) && (para->getParH(level)->coordY_SP[i] <= this->endCoordY)) && 
//			((para->getParH(level)->coordZ_SP[i] >= this->startCoordZ) && (para->getParH(level)->coordZ_SP[i] <= this->endCoordZ)) )
//		{
//			if (para->getParH(level)->geoSP[i] >= GEO_FLUID)
//			{
//				para->getParH(level)->geoSP[i] = this->geoID;
//				nodeIDsPorousMedia.push_back(i);
//				counter++;
//			}
//		}
//	}
//
//	this->sizePM = counter;
//
//	for (unsigned int j = 0; j <= this->sizePM; j++)
//	{
//	}
//}

//getter
double PorousMedia::getPorosity(){					return this->porosity; }
double PorousMedia::getDarcySI() {					return this->darcySI; }
double PorousMedia::getForchheimerSI(){				return this->forchheimerSI; }
double PorousMedia::getDarcyLBM() {					return this->darcyLBM; }
double PorousMedia::getForchheimerLBM() {			return this->forchheimerLBM; }
unsigned int PorousMedia::getGeoID() {				return this->geoID; }
unsigned int PorousMedia::getSizePM() {				return this->sizePM; }
unsigned int* PorousMedia::getHostNodeIDsPM() {		return this->hostNodeIDsPM; }
unsigned int* PorousMedia::getDeviceNodeIDsPM() {	return this->deviceNodeIDsPM; }
