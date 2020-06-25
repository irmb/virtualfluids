#pragma once
#include "Standard.h"
#include "CoordNeighborGeoV.h"

using namespace std;

class MeasuredPoints 
	: public CoordNeighborGeoV 
{
private:
	
public:
	MeasuredPoints(void);
	MeasuredPoints(string ad);
	~MeasuredPoints(void);

	void init();

};

