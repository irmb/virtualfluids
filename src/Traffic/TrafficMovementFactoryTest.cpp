#include "TrafficMovementFactoryTest.h"

#include <iostream>

#include "GridGenerator/StreetPointFinder/JunctionReader.h"
#include "GridGenerator/StreetPointFinder/SourceReader.h"
#include "GridGenerator/StreetPointFinder/SinkReader.h"

#include "RoadNetwork/RoadMaker.h"
#include "TrafficMovement.h"
#include "Source/SourceRandom.h"
#include "Junction/JunctionRandom.h"
#include "Sink/SinkRandom.h"
#include "Output/ConcentrationByPosition.h"
#include "Output/ConcBySpeedAndAcceleration.h"
#include "Utilities/safe_casting.h"


void TrafficMovementFactoryTest::initTrafficMovement(real * pconcArrayStart)
{
	//Variables

	real vehicleDensity = 0.1f;

	const uint maxVelocity = 14;
	real dawdlePossibility = (real) 0.2; //typical value: 0.2
	real slowToStartPossibility = (real) 0.4;

	uint vehicleLength = 7;

	uint roadLength = 10;


	//make RoadNetwork
	auto roadNetwork = std::make_unique<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);


	//Sources
	std::unique_ptr<Source> source = std::make_unique <SourceRandom>(SourceRandom(0, 0.7f, maxVelocity));
	roadNetwork->addSource(source);


	//Sinks
	std::unique_ptr<Sink> s = std::make_unique <SinkRandom>(SinkRandom(9, 0.5f));
	roadNetwork->addSink(move(s));


	//init TrafficMovement
	this->simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	simulator->setSlowToStart(slowToStartPossibility);
	simulator->setMaxAcceleration(2);


	////init ConcentrationOutwriter
	//std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), pconcArrayStart));
	//simulator->setConcentrationOutwriter(move(writer));
}

void TrafficMovementFactoryTest::calculateTimestep(uint step, uint stepForVTK)
{
	simulator->calculateTimestep(step);
}

void TrafficMovementFactoryTest::loopThroughTimesteps(uint timeSteps)
{
	simulator->setSaveResultsTrue(timeSteps);
	simulator->loopTroughTimesteps(timeSteps);
}
