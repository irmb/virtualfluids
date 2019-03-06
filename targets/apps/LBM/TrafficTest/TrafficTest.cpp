#include <iostream>
#include <vector>
#include <memory>

#include "Core/DataTypes.h"

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"
#include "GridGenerator/StreetPointFinder/JunctionReader.h"
#include "GridGenerator/StreetPointFinder/SourceReader.h"
#include "GridGenerator/StreetPointFinder/SinkReader.h"

#include "Traffic/RoadNetwork/RoadMaker.h"
#include "Traffic/TrafficMovement.h"
#include "Traffic/Source/SourceRandom.h"
#include "Traffic/Junction/Junction.h"
#include "Traffic/Junction/JunctionRandom.h"
#include "Traffic/Sink/Sink.h"
#include "Traffic/Sink/SinkRandom.h"
#include "Traffic/Output/ConcentrationByPosition.h"
#include "Traffic/Output/ConcBySpeedAndAcceleration.h"


int main()
{
	//Variables

	int numberOfTimesteps = 1000;
	float vehicleDensity = 0.5; 

	const uint maxVelocity = 14;
	float dawdlePossibility = (float) 0.2; //typical value: 0.2
	float slowToStartPossibility = (float) 0.4;

	uint vehicleLength = 7;


	//Logger

	logging::Logger::addStream(&std::cout);
	logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
	logging::Logger::timeStamp(logging::Logger::ENABLE);
	logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);


	//StreetPointFinder

	StreetPointFinder finder;

	finder.readStreets("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	finder.writeVTK("C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/ExampleStreets.vtk");

	JunctionReader junctionReader;
	junctionReader.readJunctions("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Junctions.txt", finder);

	SinkReader sinkReader;
	sinkReader.readSinks("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sinks.txt", finder);

	SourceReader sourceReader;
	sourceReader.readSources("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sources.txt", finder);


	{ 
		//calculate RoadLength
		uint roadLength = 0;
		uint numberOfStreets = finder.streets.size();
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
		auto simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
		//simulator->setSaveResultsTrue();
		simulator->setSlowToStart(slowToStartPossibility);
		simulator->initCalculation(numberOfTimesteps);


		//init ConcentrationOutwriter
		//std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcentrationByPosition>(ConcentrationByPosition(simulator->getRoadLength()));
		//simulator->setConcentrationOutwriter(move(writer));


		//prepare writing to vtk
		std::string outputPath("C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/");
		std::string outputFilename("ExampleStreets");
		auto& cars = simulator->getVehiclesForVTK();


		//loop through timesteps
		for (int step = 1; step < numberOfTimesteps + 1; step++) {
			simulator->calculateTimestep(step);
			simulator->visualizeVehicleLengthForVTK();
			finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(step) + ".vtk", cars);
		}

	
		std::vector<float> sim = simulator->getConcentrations();
		simulator->dispResults();

		std::cout << std::endl << std::endl;
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//std::vector<int> initialDistribution = { 1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1 };
	//std::vector<int> oneCar = { -1,1,-1,-1,-1,-1 };
	//std::vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, 2, -1, -1,-1,0,-1,-1,-1,-1,1,-1,-1,-1 };
	//int roadLength = 20;
	//vehicleLength = 3;


	////OneWay random  concwriter-test
	//{
	//	std::cout << "OneWay random" << std::endl;


	//	auto roadNetwork    = std::make_unique<RoadMak er>(oneCar, maxVelocity, vehicleLength);
	//	auto simulator = std::make_shared<TrafficMovement>(move(roadNetwork) ,dawdlePossibility);

	//	auto writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), 5, 0));
	//	simulator->setSaveResultsTrue();
	//	simulator->setConcentrationOutwriter(move(writer));

	//	simulator->initCalculation(1);
	//	simulator->loopTroughTimesteps();

	//	std::cout << std::endl << std::endl;

	//}


	////JunctionTest
	//{
	//	auto roadNetwork = std::make_unique<RoadMaker>(fiveCars, maxVelocity, vehicleLength);

	//	vector <uint> in4 = { 9 };
	//	vector<uint> out4 = { 10 };
	//	std::unique_ptr<Junction>  j = std::make_unique<JunctionRandom>(in4, out4);
	//	j->setCellIndexForNoUTurn({ 10 });
	//	roadNetwork->addJunction(j);

	//	std::shared_ptr<TrafficMovement> simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	//	simulator->setSaveResultsTrue();
	//	simulator->initCalculation(numberOfTimesteps);

	//	for (int step = 1; step < numberOfTimesteps + 1; step++) {
	//		simulator->calculateTimestep(step);
	//	}

	//	std::cout << "Number of Cars: " << simulator->getNumberOfCars() << std::endl;
	//	simulator->dispResults();
	//}

	//{
	//	dawdlePossibility = 0;
	//	vector<int> car = { -1,-1,-1,4,-1,-1,-1,-1,-1 };

	//	unique_ptr<RoadMaker> r = make_unique<RoadMaker>(car, 5, 2);

	//	vector <uint> in = { 4 };
	//	vector <uint> out = { 5 };
	//	unique_ptr<Junction>  j = make_unique<JunctionRandom>(JunctionRandom(in, out));
	//	r->addJunction(j);

	//	shared_ptr<TrafficMovement> simulator = make_shared<TrafficMovement>(move(r), dawdlePossibility);
	//	simulator->setSaveResultsTrue();
	//	simulator->initCalculation(1);

	//	simulator->calculateTimestep(1);

	//	bool success = false;
	//	if (simulator->getSpeedAtPosition(7) == 5) success = true;
	//	std::cout << success << std::endl;
	//	simulator->dispResults();
	//}



	////OneWay initial distribution
	//{
	//	cout << "OneWay initial distribution" << endl;
	//	unique_ptr<RoadNetworkData> roadNetwork = make_unique<RoadMaker>(initialDistribution, maxVelocity, vehicleLength);
	//	shared_ptr<TrafficMovement> simulator = make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	//	simulator->initCalculation(numberOfTimesteps);
	//simulator->loopTroughTimesteps();
	//	simulator->dispResults();	
	//	cout << endl << endl;
	//}


	//////OneWay slowToStart
	//{
	//	dawdlePossibility = 0;
	//	cout << "OneWay slowToStart" << endl;
	//	vector<int> initialDist = { -1,1,5,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,0,-1,-1,-1,-1,-1 };

	//	unique_ptr<RoadNetworkData> roadNetwork = make_unique<RoadMaker>(initialDist, maxVelocity, 1);
	//	shared_ptr<TrafficMovement> simulator = make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);
	//	simulator->setSlowToStart(static_cast<float>(0.9999));
	//	simulator->initCalculation(numberOfTimesteps);
	//simulator->loopTroughTimesteps();
	//	simulator->dispResults();
	//	cout << endl << endl;
	//}



	////sources and sinks
	//{
	//	cout << "sources and sinks" << endl;

	//	vector< unique_ptr<Sink> > sinks;
	//	sinks.push_back(make_unique <SinkRandom>(5, static_cast<float>(0.5)));

	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(oneCar, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceRandom>(0, static_cast<float>(1.0), roadNetwork->getMaxVelocity()));

	//	roadNetwork->setSources(move(sources));
	//	roadNetwork->setSinks(move(sinks));

	//	shared_ptr<TrafficMovement> roadSource = make_shared<TrafficMovement>(move(roadNetwork), static_cast<float>(0.0));

	//	roadSource->initCalculation(numberOfTimesteps);
	//	roadSource->loopTroughTimesteps();
	//	roadSource->dispResults();
	//	cout << endl << endl;
	//}



	////mergingRoad
	//{
	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(25, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceRandom>(0, static_cast<float>(1), roadNetwork->getMaxVelocity()));
	//	sources.push_back(make_unique <SourceRandom>(10, static_cast<float>(1), roadNetwork->getMaxVelocity()));
	//	roadNetwork->setSources(move(sources));

	//	unique_ptr<Sink> s = make_unique <SinkRandom>(SinkRandom(24, static_cast<float>(0.0)));
	//	roadNetwork->addSink(move(s));

	//	vector< unique_ptr<Junction> > junctions;
	//	vector<uint> in = { 9,20 };
	//	vector<uint> out = { 21 };
	//	junctions.push_back(make_unique <JunctionRandom>(JunctionRandom(in, out)));

	//	roadNetwork->setJunctions(junctions);

	//	shared_ptr<TrafficMovement> mergingRoad = make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);


	//	mergingRoad->initCalculation(numberOfTimesteps);
	//	mergingRoad->loopTroughTimesteps();
	//	mergingRoad->dispResults();
	//	cout << endl << endl;
	//}
}

