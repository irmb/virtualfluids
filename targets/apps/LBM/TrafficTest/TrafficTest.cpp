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


static std::shared_ptr<TrafficMovement> initTrafficSim(StreetPointFinder &finder, float * pconcArrayStart = nullptr, uint concArraySize = 0) {

	//Variables

	float vehicleDensity = 0.1;

	const uint maxVelocity = 14;
	float dawdlePossibility = (float) 0.2; //typical value: 0.2
	float slowToStartPossibility = (float) 0.4;

	uint vehicleLength = 7;


	//StreetPointFinder

	finder.readStreets("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	finder.writeVTK("C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/ExampleStreets.vtk");

	JunctionReader junctionReader;
	junctionReader.readJunctions("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Junctions.txt", finder);

	SinkReader sinkReader;
	sinkReader.readSinks("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sinks.txt", finder);

	SourceReader sourceReader;
	sourceReader.readSources("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/Basel/resources/Sources.txt", finder);

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
	simulator = simulator;
	simulator->setSlowToStart(slowToStartPossibility);


	//init ConcentrationOutwriter
	std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), simulator->getMaxVelocity()));
	//std::unique_ptr<ConcentrationOutwriter> writer = std::make_unique<ConcentrationByPosition>(ConcentrationByPosition(simulator->getRoadLength()), pconcArrayStart, concArraySize);
	simulator->setConcentrationOutwriter(move(writer));
	return simulator;

}


static void calculateTimestep(std::shared_ptr<TrafficMovement>& simulator, uint step) {
	simulator->calculateTimestep(step);
};

int main()
{
	{
		uint numberOfTimesteps = 500;


		//Logger

		logging::Logger::addStream(&std::cout);
		logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
		logging::Logger::timeStamp(logging::Logger::ENABLE);
		logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);

		StreetPointFinder finder;
		std::shared_ptr<TrafficMovement> simulator = initTrafficSim(finder);


		////save results
		//simulator->setSaveResultsTrue(numberOfTimesteps);


		//prepare writing to vtk
		std::string outputPath("C:/Users/hiwi/BaselDokumente/Basel_Ergebnisse/");
		std::string outputFilename("ExampleStreets");
		auto& cars = simulator->getVehiclesForVTK();


		//loop through timesteps

		for (uint step = 1; step < numberOfTimesteps + 1; step++) {
			calculateTimestep(simulator, step);
			simulator->visualizeVehicleLengthForVTK();
			finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(step) + ".vtk", cars);
		}

		//simulator->dispResults();

		std::cout << std::endl << std::endl;

	}








	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	//std::vector<int> oneCar = { -1,1,-1,-1,-1,-1 };
	//std::vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, -1, 1, -1,-1,-1,-1,0,-1,-1,1,-1,-1,-1 };
	//int roadLength = 20;
	//uint vehicleLength = 2;
	//uint maxVelocity = 5;
	//float dawdlePossibility = 0;

	//uint numberOfTimesteps = 10;

	////concwriter-test
	//{
	//	
	//	std::vector<int> cars = { 1,-1,0,-1,0,-1,-1,-1,-1,-1 }; //IdleTest
	//	//cars = { -1,-1,0,-1,-1,1,-1,0,-1,0 }; //IdleTest2

	//	std::cout << "concwriter-test" << std::endl;


	//	auto roadNetwork = std::make_unique<RoadMaker>(cars, maxVelocity, vehicleLength);
	//	auto simulator = std::make_shared<TrafficMovement>(move(roadNetwork), dawdlePossibility);


	//	float concArray[10];
	//	float * pconcArrayStart = &concArray[0];

	//	auto writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), pconcArrayStart, 10, maxVelocity));
	//	simulator->setSaveResultsTrue(numberOfTimesteps);
	//	simulator->setConcentrationOutwriter(move(writer));

	//	simulator->loopTroughTimesteps(numberOfTimesteps);

	//	std::cout << std::endl << std::endl;
	//}


	////JunctionTest
	//{
	//	std::cout << "junction-test" << std::endl;

	//	auto roadNetwork = std::make_unique<RoadMaker>(fiveCars, maxVelocity, vehicleLength);

	//	std::vector <uint> in4 = { 9 };
	//	std::vector<uint> out4 = { 10 };
	//	std::unique_ptr<Junction>  j = std::make_unique<JunctionRandom>(in4, out4);
	//	roadNetwork->addJunction(j);

	//	std::shared_ptr<TrafficMovement> simulator = std::make_shared<TrafficMovement>(move(roadNetwork), 0);
	//	auto writer = std::make_unique<ConcBySpeedAndAcceleration>(ConcBySpeedAndAcceleration(simulator->getRoadLength(), maxVelocity));
	//	simulator->setSaveResultsTrue(numberOfTimesteps);
	//	simulator->setConcentrationOutwriter(move(writer));

	//	simulator->loopTroughTimesteps(numberOfTimesteps);

	//	std::cout << "Number of Cars: " << simulator->getNumberOfCars() << std::endl;
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
	//	simulator->setSlowToStart(0.9999f);
	//	simulator->initCalculation(numberOfTimesteps);
	//simulator->loopTroughTimesteps();
	//	simulator->dispResults();
	//	cout << endl << endl;
	//}



	////sources and sinks
	//{
	//	cout << "sources and sinks" << endl;

	//	vector< unique_ptr<Sink> > sinks;
	//	sinks.push_back(make_unique <SinkRandom>(5, 0.5f));

	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(oneCar, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceRandom>(0,1.0f, roadNetwork->getMaxVelocity()));

	//	roadNetwork->setSources(move(sources));
	//	roadNetwork->setSinks(move(sinks));

	//	shared_ptr<TrafficMovement> roadSource = make_shared<TrafficMovement>(move(roadNetwork), 0.0f);

	//	roadSource->initCalculation(numberOfTimesteps);
	//	roadSource->loopTroughTimesteps();
	//	roadSource->dispResults();
	//	cout << endl << endl;
	//}



	////mergingRoad
	//{
	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(25, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceRandom>(0, 1f, roadNetwork->getMaxVelocity()));
	//	sources.push_back(make_unique <SourceRandom>(10, 1f, roadNetwork->getMaxVelocity()));
	//	roadNetwork->setSources(move(sources));

	//	unique_ptr<Sink> s = make_unique <SinkRandom>(SinkRandom(24, 0f));
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

