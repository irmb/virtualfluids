#ifndef SimulationFileNames_H
#define SimulationFileNames_H

#include <string>
#include "GridGenerator_EXPORT.h"

struct GridGenerator_EXPORT simulationFileNames
{
	static const std::string coordX;
	static const std::string coordY;
	static const std::string coordZ;
	static const std::string neighborX;
	static const std::string neighborY;
	static const std::string neighborZ;
	static const std::string geoVec;
	
	static const std::string geomBoundaryQ;
	static const std::string geomBoundaryValues;
	
	static const std::string topBoundaryQ;
	static const std::string topBoundaryValues;
	
	static const std::string bottomBoundaryQ;
	static const std::string bottomBoundaryValues;
	
	static const std::string frontBoundaryQ;
	static const std::string frontBoundaryValues;
	
	static const std::string backBoundaryQ;
	static const std::string backBoundaryValues;
	
	static const std::string inletBoundaryQ;
	static const std::string inletBoundaryValues;
	
	static const std::string outletBoundaryQ;
	static const std::string outletBoundaryValues;
};


#endif
