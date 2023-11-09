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

PorousMedia::PorousMedia(double porosity, unsigned int geoID, double darcySI, double forchheimerSI, double dxLBM, double dtLBM, unsigned int level)
{
	this->porosity = porosity;
	this->geoID = geoID;
	this->darcySI = darcySI;
	this->forchheimerSI = forchheimerSI;
	this->dxLBM = dxLBM;
	this->dtLBM = dtLBM;
	//this->lengthOfPorousMedia = 1.0;
	//this->forchheimerLBM = this->forchheimerSI * this->dxLBM / this->lengthOfPorousMedia;
	//this->darcyLBM = this->darcySI * this->dtLBM / this->lengthOfPorousMedia;
	this->level = level;
	setCoordinatesToZero();
	setResistanceLBM();
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

void PorousMedia::setSizePM(unsigned int sizePM)
{
	this->sizePM = sizePM;
}

void PorousMedia::setResistanceLBM() 
{
	if ((this->endCoordX - this->startCoordX) > 1.0)
	{
		this->lengthOfPorousMedia = (this->endCoordX - this->startCoordX);
	}
	else
	{
		this->lengthOfPorousMedia = 1.0;
	}
	this->forchheimerLBM = this->forchheimerSI * this->dxLBM / this->lengthOfPorousMedia;
	this->darcyLBM = this->darcySI * this->dtLBM / this->lengthOfPorousMedia;
}

//void PorousMedia::definePMarea(Parameter* para, unsigned int level)
//{
//	unsigned int counter = 0;
//	for (unsigned int i = 0; i < para->getParH(level)->numberOfNodes; i++)
//	{
//		if (((para->getParH(level)->coordinateX[i] >= this->startCoordX) && (para->getParH(level)->coordinateX[i] <= this->endCoordX)) &&
//			((para->getParH(level)->coordinateY[i] >= this->startCoordY) && (para->getParH(level)->coordinateY[i] <= this->endCoordY)) && 
//			((para->getParH(level)->coordinateZ[i] >= this->startCoordZ) && (para->getParH(level)->coordinateZ[i] <= this->endCoordZ)) )
//		{
//			if (para->getParH(level)->typeOfGridNode[i] >= GEO_FLUID)
//			{
//				para->getParH(level)->typeOfGridNode[i] = this->geoID;
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
unsigned int PorousMedia::getLevelPM() {			return this->level; }
unsigned int* PorousMedia::getHostNodeIDsPM() {		return this->hostNodeIDsPM; }
unsigned int* PorousMedia::getDeviceNodeIDsPM() {	return this->deviceNodeIDsPM; }
double PorousMedia::getStartX(){					return this->startCoordX; }
double PorousMedia::getStartY(){					return this->startCoordY; }
double PorousMedia::getStartZ(){					return this->startCoordZ; }
double PorousMedia::getEndX(){						return this->endCoordX; }
double PorousMedia::getEndY(){						return this->endCoordY; }
double PorousMedia::getEndZ(){						return this->endCoordZ; }

