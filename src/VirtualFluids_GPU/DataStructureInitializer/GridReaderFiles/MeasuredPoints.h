#pragma once

#include <string>
#include <vector>

#include "CoordNeighborGeoV.h"


class MeasuredPoints 
	: public CoordNeighborGeoV 
{
private:
    std::vector< std::vector<unsigned int> > points;

public:
	MeasuredPoints();
	MeasuredPoints(std::string ad);
	~MeasuredPoints();

	void init();

};

