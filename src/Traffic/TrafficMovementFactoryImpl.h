# pragma once

//#include <VirtualFluidsDefinitions.h>
#include "TrafficMovementFactory.h"

#include <vector>
#include <memory>

#include "Core/DataTypes.h"
#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

class TrafficMovement;

class VF_PUBLIC TrafficMovementFactoryImpl:
public TrafficMovementFactory{
public:
	TrafficMovementFactoryImpl();
	~TrafficMovementFactoryImpl() {};
	virtual void initTrafficMovement(real * pconcArrayStart = nullptr);
	virtual void calculateTimestep(uint step, uint stepForVTK);

protected:
	StreetPointFinder finder;
	std::shared_ptr<TrafficMovement> simulator;
	std::string outputPath;
	std::string outputFilename;
	const std::vector<int>* cars;

};