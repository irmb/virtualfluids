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

	getRoadLength = 10;


	//make RoadNetwork
	auto roadNetwork = std::make_unique<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);


	//Sources
	unique_ptr<Source> s = make_unique <SourceRandom>(SourceRandom(9, 0.5f))
	roadNetwork->addSource(source);


	//Sinks
	unique_ptr<Sink> s = make_unique <SinkRandom>(SinkRandom(9, 0.5f));
	/ roadNetwork->addSink(move(s));


	//init TrafficMovement
	this->simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	simulator->setSlowToStart(slowToStartPossibility);
	simulator->setMaxAcceleration(2);
	simulator->se


	//init ConcentrationOutwriter
	std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), pconcArrayStart));
	simulator->setConcentrationOutwriter(move(writer));

	//prepare writing to vtk
	//this->outputPath = "M:/Basel2019/results/";
	this->outputPath = "C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/";
	this->outputFilename = "ExampleStreets";
	this->cars = &(simulator->getVehiclesForVTK());

	//write initial Timestep
	simulator->visualizeVehicleLengthForVTK();
	finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(0) + ".vtk", *cars);
}

void TrafficMovementFactoryTest::calculateTimestep(uint step, uint stepForVTK)
{
	simulator->calculateTimestep(step);
	simulator->visualizeVehicleLengthForVTK();
	finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(stepForVTK) + ".vtk", *cars);
}

void TrafficMovementFactoryTest::loopThroughTimesteps(uint timeSteps)
{
	simulator->setSaveResultsTrue(timeSteps);
	calculateTimeStep;
}
