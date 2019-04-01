# pragma once

#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

#include "Core/DataTypes.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

class TrafficMovement;

class VF_PUBLIC TrafficMovementFactory
{
public:
	TrafficMovementFactory();
	~TrafficMovementFactory() {};
	virtual void initTrafficMovement(std::string path, real * pConcArray = nullptr);
	virtual void calculateTimestep(uint step, uint stepForVTK);

protected:
	StreetPointFinder finder;
	std::shared_ptr<TrafficMovement> simulator;
	std::string inputPath;
	std::string outputPath;
	std::string outputFilename;
	const std::vector<int>* cars;

};