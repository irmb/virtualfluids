# pragma once

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include <vector>
#include <memory>

#include "Traffic_export.h"

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

class TrafficMovement;
class TrafficLogger;

class TRAFFIC_EXPORT TrafficMovementFactory
{
public:
	TrafficMovementFactory();
	~TrafficMovementFactory() {};
	virtual void initTrafficMovement(std::string path, bool useGPU, real * pConcArray = nullptr, int* naschVelocity = nullptr);
	virtual void calculateTimestep(uint step);
	virtual void writeTimestep(uint stepForVTK);
	void writeReducedTimestep(uint stepForVTK);
	virtual void endSimulation(uint numTimesteps, double duration);

protected:
	StreetPointFinder finder;
	std::shared_ptr<TrafficMovement> simulator;

	std::string inputPath;
	std::string outputPath;
	std::string outputFilename;
	const std::vector<int>* cars;

	bool useLogger;
	bool useGPU;
};