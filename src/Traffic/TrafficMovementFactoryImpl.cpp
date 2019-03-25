#include "TrafficMovementFactoryImpl.h"

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


TrafficMovementFactoryImpl::TrafficMovementFactoryImpl()
{
}

void TrafficMovementFactoryImpl::initTrafficMovement(real * pconcArrayStart)
{
	//Variables

	real vehicleDensity = 0.1f;

	const uint maxVelocity = 14;
	real dawdlePossibility = (real) 0.2; //typical value: 0.2
	real slowToStartPossibility = (real) 0.4;

	uint vehicleLength = 7;


	//StreetPointFinder M:\Basel2019  C:\Users\schoen\Desktop\git\MS2
	//finder.readStreets("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	//finder.writeVTK("M:/Basel2019/results/ExampleStreets.vtk");
	finder.readStreets("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	finder.writeVTK("C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/ExampleStreets.vtk");

	JunctionReader junctionReader;
	//junctionReader.readJunctions("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Junctions.txt", finder);
	junctionReader.readJunctions("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Junctions.txt", finder);

	SinkReader sinkReader;
	//sinkReader.readSinks("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Sinks.txt", finder);
	sinkReader.readSinks("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sinks.txt", finder);

	SourceReader sourceReader;
	//sourceReader.readSources("C:/Users/schoen/Desktop/git/MS2/git/targets/apps/LBM/Basel/resources/Sources.txt", finder);
	sourceReader.readSources("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sources.txt", finder);

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
		junctions.push_back(std::make_unique <JunctionRandom>(junctionReader.junctions[i].inCells, junctionReader.junctions[i].outCells));
		junctions[i]->setCellIndexForNoUTurn(junctionReader.junctions[i].carCanNotEnterThisOutCell);
	}
	roadNetwork->setJunctions(move(junctions));


	//init TrafficMovement
	this->simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	simulator->setSlowToStart(slowToStartPossibility);
	simulator->setMaxAcceleration(2);
	simulator->setUseGPU();

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

void TrafficMovementFactoryImpl::calculateTimestep(uint step, uint stepForVTK)
{
	simulator->calculateTimestep(step);
	simulator->visualizeVehicleLengthForVTK();
	finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(stepForVTK) + ".vtk", *cars);
}

