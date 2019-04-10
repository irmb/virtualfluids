#include "TrafficMovementFactory.h"

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
#include "Utilities/Logger.h"


TrafficMovementFactory::TrafficMovementFactory()
{
}

void TrafficMovementFactory::initTrafficMovement(std::string path, real * pConcArray)
{
	//Variables

	real vehicleDensity = 0.01f;

	uint vehicleLength = 7;
	uint maxVelocity = 14;
	uint maxAcceleration = 2;

	real dawdlePossibility = (real) 0.2; //typical value: 0.2
	real slowToStartPossibility = (real) 0.3;

	bool useGPU = true;
	bool useSlowToStart = true;
	useLogger = true;

	std::string info = "Only Traffic, writing vtk";
	


	//Paths

	inputPath = path + "VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/";
	outputPath = path + "Basel_Ergebnisse/";
	outputFilename = "Basel_Traffic";
	std::string logfile = outputPath + "TrafficLog.txt";



	//TrafficLogger	
	if (useLogger) {
		TrafficLogger::startLogger(logfile);
		TrafficLogger::writeSimulationStart(info, useGPU);
	}


	//StreetPointFinder M:\Basel2019  C:\Users\schoen\Desktop\git\MS2
	//finder.readStreets("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	//finder.writeVTK("M:/Basel2019/results/ExampleStreets.vtk");
	finder.readStreets(inputPath + "Streets.txt");
	finder.writeVTK(outputPath + outputFilename + ".vtk");


	JunctionReader junctionReader;
	//junctionReader.readJunctions("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Junctions.txt", finder);
	junctionReader.readJunctions(inputPath + "Junctions.txt", finder);


	SinkReader sinkReader;
	//sinkReader.readSinks("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Sinks.txt", finder);
	sinkReader.readSinks(inputPath + "Sinks.txt", finder);


	SourceReader sourceReader;
	//sourceReader.readSources("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Sources.txt", finder);
	sourceReader.readSources(inputPath + "Sources.txt", finder);


	//calculate RoadLength
	uint roadLength = 0;
	uint numberOfStreets = castSizeT_Uint(finder.streets.size());
	for (uint i = 0; i < numberOfStreets; i++) {
		roadLength += finder.streets[i].numberOfCells;
	}


	//make RoadNetwork
	auto roadNetwork = std::make_unique<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);


	//Sources
	std::vector< std::unique_ptr<Source> > sources;
	for (uint i = 0; i < sourceReader.sources.size(); i++)
		sources.push_back(std::make_unique <SourceRandom>(sourceReader.sources[i].sourceIndex, sourceReader.sources[i].sourcePossibility, roadNetwork->getMaxVelocity()));
	roadNetwork->setSources(move(sources));


	//Sinks
	std::vector< std::unique_ptr<Sink> > sinks;
	for (uint i = 0; i < sinkReader.sinks.size(); i++)
		sinks.push_back(std::make_unique <SinkRandom>(sinkReader.sinks[i].sinkIndex, sinkReader.sinks[i].sinkBlockedPossibility));
	roadNetwork->setSinks(move(sinks));


	//Junctions
	std::vector <std::unique_ptr<Junction> > junctions;
	for (uint i = 0; i < junctionReader.junctions.size(); i++) {
		junctions.push_back(std::make_unique <JunctionRandom>(junctionReader.junctions[i].inCells, junctionReader.junctions[i].outCells, junctionReader.junctions[i].trafficLightSwitchTime));
		junctions[i]->setCellIndexForNoUTurn(junctionReader.junctions[i].carCanNotEnterThisOutCell);
	}
	roadNetwork->setJunctions(move(junctions));


	//init TrafficMovement
	this->simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	simulator->setMaxAcceleration(maxAcceleration);
	if (useSlowToStart) simulator->setSlowToStart(slowToStartPossibility);
	if (useLogger) simulator->setUseLogger();


	//init ConcentrationOutwriter
	if (!useGPU) {
		std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), pConcArray));
		simulator->setConcentrationOutwriter(move(writer));
	}


	//prepare writing to vtk
	//this->outputPath = "M:/Basel2019/results/";
	this->cars = &(simulator->getVehiclesForVTK());


	//write initial Timestep
	simulator->visualizeVehicleLengthForVTK();
	finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(0) + ".vtk", *cars);


	//GPU
	if (useGPU) simulator->setUseGPU(pConcArray);
}


void TrafficMovementFactory::calculateTimestep(uint step)
{
	simulator->calculateTimestep(step);

	//std::cout << "Number of Cars: " << simulator->getNumberOfCars() << std::endl;
}

void TrafficMovementFactory::writeTimestep(uint stepForVTK)
{
	simulator->visualizeVehicleLengthForVTK();
	finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(stepForVTK) + ".vtk", *cars);
}

void TrafficMovementFactory::endSimulation(uint numTimesteps, double duration)
{
	if (!useLogger) return;
	TrafficLogger::writeSimulationEnd(simulator->getRoadLength(), numTimesteps, duration);
}



