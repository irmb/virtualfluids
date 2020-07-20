#pragma once

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include <string>
#include <fstream>

class VIRTUALFLUIDS_GPU_EXPORT TrafficLogger
{
private:
	std::string filename;	
	std::ofstream file;
	static TrafficLogger instance;

public:	
	TrafficLogger() {};
	TrafficLogger(const TrafficLogger& logger) {}
	static void startLogger(std::string filename);

	static void writeSimulationStart(std::string info, bool useGPU);
	static void writeError(std::string error, uint currentTimestep);
	static void writeSimulationEnd(uint numRoadCells, uint numTimesteps, double duration);	
};

