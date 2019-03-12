# pragma once
//#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <memory>

#include "Core/DataTypes.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

class TrafficMovement;
//struct StreetPointFinder;

class VF_PUBLIC TrafficMovementFactory {
public:
	TrafficMovementFactory();
	~TrafficMovementFactory() {};
	void initTrafficMovement(real * pconcArrayStart = nullptr);
	void calculateTimestep(uint step, uint stepForVTK);


private:
	StreetPointFinder finder;
	std::shared_ptr<TrafficMovement> simulator;
	std::string outputPath;
	std::string outputFilename;
	const std::vector<int>* cars;

};