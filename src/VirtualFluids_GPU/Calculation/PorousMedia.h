#ifndef POROUS_MEDIA_H
#define POROUS_MEDIA_H

//#include "Parameter/Parameter.h"
//#include "Utilities/StringUtil.hpp"
//#include "basics/utilities/UbSystem.h"

#include <iostream>
#include <vector>
#include <stdio.h>
using namespace std;

class PorousMedia
{
public:
	PorousMedia();
	PorousMedia(double porosity, unsigned int geoID, double darcySI, double forchheimerSI, double dxLBM, double dtLBM);
	~PorousMedia();

	//setter
	void setPorosity(double porosity);
	void setDarcySI(double darcySI);
	void setForchheimerSI(double forchheimerSI);
	void setStartCoordinates(double startCoordX, double startCoordY, double startCoordZ);
	void setEndCoordinates(double endCoordX, double endCoordY, double endCoordZ);
	//void definePMarea(Parameter* para, unsigned int level);
	void setHostNodeIDsPM(unsigned int* hostNodeIDsPM);
	void setDeviceNodeIDsPM(unsigned int* deviceNodeIDsPM);

	//getter
	double getPorosity();
	double getDarcySI();
	double getForchheimerSI();
	double getDarcyLBM();
	double getForchheimerLBM();
	unsigned int getGeoID();
	unsigned int getSizePM();
	unsigned int* getHostNodeIDsPM();
	unsigned int* getDeviceNodeIDsPM();

private:
	double porosity;
	double darcySI; //[1/s]
	double darcyLBM; 
	double forchheimerSI; //[1/m]
	double forchheimerLBM;
	double dxLBM;
	double dtLBM;
	double startCoordX;
	double startCoordY;
	double startCoordZ;
	double endCoordX;
	double endCoordY;
	double endCoordZ;
	unsigned int geoID;
	std::vector< unsigned int > nodeIDsPorousMedia;
	unsigned int sizePM;
	unsigned int *hostNodeIDsPM;
	unsigned int *deviceNodeIDsPM;

	void setCoordinatesToZero() 
	{
		startCoordX = 0;
		startCoordY = 0;
		startCoordZ = 0;
		endCoordX = 0;
		endCoordY = 0;
		endCoordZ = 0;
	};

	void setSIvaluesToZero()
	{
		darcySI = 0.0;
		forchheimerSI = 0.0;
	};

	void setLBMvaluesToZero()
	{
		darcyLBM = 0.0;
		forchheimerLBM = 0.0;
		dxLBM = 0.0;
		dtLBM = 0.0;
	};


};


#endif /* POROUS_MEDIA_H */
