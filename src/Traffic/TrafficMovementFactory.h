# pragma once
//#include <VirtualFluidsDefinitions.h>

#include <iostream>
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
	void initTrafficMovement(real * pconcArrayStart = nullptr, uint concArraySize = 0);
	void calculateTimestep(uint step);


private:
	StreetPointFinder finder;
	std::shared_ptr<TrafficMovement> simulator;
	std::string outputPath;
	std::string outputFilename;
	const std::vector<int>* cars;

};