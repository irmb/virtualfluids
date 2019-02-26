#include <iostream>
#include <vector>
#include <memory>

#include "GridGenerator/StreetPointFinder/StreetPointFinder.h"

#include "Traffic/RoadMaker.h"
#include "Traffic/OneWayRoad.h"
#include "Traffic/OneWayRoadSSJ.h"
#include "Traffic/SourceTerm.h"
#include "Traffic/Junction.h"
#include "Traffic/JunctionSimple.h"
#include "Traffic/Sink.h"
#include "Traffic/SinkSimple.h"
#include "Traffic/ConcentrationByPosition.h"

//using namespace std;


int main()
{
	//Logger

	logging::Logger::addStream(&std::cout);
	logging::Logger::setDebugLevel(logging::Logger::Level::INFO_LOW);
	logging::Logger::timeStamp(logging::Logger::ENABLE);
	logging::Logger::enablePrintedRankNumbers(logging::Logger::ENABLE);




	//Variables

	int numberOfTimesteps = 100;
	float vehicleDensity = 0.05;

	const unsigned int maxVelocity = 5;
	float dawdlePossibility = 0.5;//typical value: 0.2
	unsigned int vehicleLength = 7;




	//StreetPointFinder

	StreetPointFinder finder;

	finder.readStreets("C:/Users/hiwi/BaselDokumente/VirtualFluidsGPU/git/targets/apps/LBM/streetTest/resources/ExampleStreets.txt");
	finder.writeVTK("C:/Users/hiwi/Desktop/Basel_Ergebnisse/ExampleStreets.vtk");





	//One RandomRoad 

	{
		std::cout << "OneWay random" << std::endl;



		unsigned int roadLength = 0;
		for (unsigned int i = 0; i < 2; i++) {
			roadLength += finder.streets[i].numberOfCells;
		}


		auto roadNetwork = std::make_unique<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);
		auto simulator = std::make_shared<OneWayRoad>(move(roadNetwork), dawdlePossibility);
		simulator->initCalculation(numberOfTimesteps);

		//unique_ptr<ConcentrationOutwriter> writer = make_unique<ConcentrationByPosition>(ConcentrationByPosition(simulator->getRoadLength()));
		//simulator->setConcentrationOutwriter(move(writer));




		std::string outputPath("C:/Users/hiwi/Desktop/Basel_Ergebnisse/");
		std::string outputFilename("ExampleStreets");

		auto& cars = simulator->getVehiclesForVTK();

		for (int step = 1; step < numberOfTimesteps + 1; step++) {
			simulator->calculateTimestep(step);
			simulator->visualizeVehicleLengthForVTK();
			finder.writeVTK(outputPath + outputFilename + "_" + std::to_string(step) + ".vtk", cars);
		}
	

		std::cout << std::endl << std::endl;

	}







	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> initialDistribution = { 1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1 };
	std::vector<int> oneCar = { 1,-1,-1,-1,-1,-1 };
	std::vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, 2, -1, -1,-1,0,-1,-1,-1,-1,1,-1,-1,-1 };
	int roadLength = 20;
	vehicleLength = 3;



	////OneWay random 
	//{
	//	std::cout << "OneWay random" << std::endl;


	//	auto roadNetwork    = std::make_unique<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);
	//	auto simulator = std::make_shared<OneWayRoad>(move(roadNetwork) ,dawdlePossibility);

	//	//auto writer = make_unique<ConcentrationByPosition>(ConcentrationByPosition(simulator->getRoadLength()));
	//	//simulator->setConcentrationOutwriter(writer);

	//	simulator->initCalculation(numberOfTimesteps);
	//simulator->loopTroughTimesteps();
	//	simulator->dispResults();

	//	std::cout << std::endl << std::endl;

	//}


	////OneWay initial distribution
	//{
	//	cout << "OneWay initial distribution" << endl;
	//	unique_ptr<RoadNetworkData> roadNetwork = make_unique<RoadMaker>(initialDistribution, maxVelocity, vehicleLength);
	//	shared_ptr<OneWayRoad> simulator = make_shared<OneWayRoad>(move(roadNetwork), dawdlePossibility);
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
	//	shared_ptr<OneWayRoad> simulator = make_shared<OneWayRoad>(move(roadNetwork), dawdlePossibility);
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
	//	sinks.push_back(make_unique <SinkSimple>(5, static_cast<float>(0.5)));

	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(oneCar, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceTerm>(0, static_cast<float>(1.0), roadNetwork->getMaxVelocity()));

	//	roadNetwork->setSources(move(sources));
	//	roadNetwork->setSinks(move(sinks));

	//	shared_ptr<OneWayRoadSSJ> roadSource = make_shared<OneWayRoadSSJ>(move(roadNetwork), static_cast<float>(0.0));

	//	roadSource->initCalculation(numberOfTimesteps);
	//	roadSource->loopTroughTimesteps();
	//	roadSource->dispResults();
	//	cout << endl << endl;
	//}



	////mergingRoad
	//{
	//	unique_ptr<RoadMaker> roadNetwork = make_unique<RoadMaker>(25, maxVelocity, vehicleLength);

	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceTerm>(0, static_cast<float>(1), roadNetwork->getMaxVelocity()));
	//	sources.push_back(make_unique <SourceTerm>(10, static_cast<float>(1), roadNetwork->getMaxVelocity()));
	//	roadNetwork->setSources(move(sources));

	//	unique_ptr<Sink> s = make_unique <SinkSimple>(SinkSimple(24, static_cast<float>(0.0)));
	//	roadNetwork->addSink(move(s));

	//	vector< unique_ptr<Junction> > junctions;
	//	vector<unsigned int> in = { 9,20 };
	//	vector<unsigned int> out = { 21 };
	//	junctions.push_back(make_unique <JunctionSimple>(JunctionSimple(in, out)));

	//	roadNetwork->setJunctions(junctions);

	//	shared_ptr<OneWayRoadSSJ> mergingRoad = make_shared<OneWayRoadSSJ>(move(roadNetwork), dawdlePossibility);


	//	mergingRoad->initCalculation(numberOfTimesteps);
	//	mergingRoad->loopTroughTimesteps();
	//	mergingRoad->dispResults();
	//	cout << endl << endl;
	//}


	////splittingRoad
	//{
	//	OneWayRoadSSJ* splittingRoad = new OneWayRoadSSJ(25, maxVelocity, dawdlePossibility, vehicleLength);

	//	vector<Source*> sources2 = { new SourceTerm(0, 1, splittingRoad) };

	//	splittingRoad->setSources(sources2);

	//	vector<Sink*> sinks2;
	//	sinks2.push_back(new SinkSimple(24, 0));
	//	sinks2.push_back(new SinkSimple(19, 0));
	//	splittingRoad->setSinks(sinks2);

	//	vector <unsigned int> in = { 9 };
	//	vector<unsigned int> out = { 11,20 };
	//	vector<Junction*> junctions2 = { new  JunctionSimple(in, out, splittingRoad) };
	//	splittingRoad->setJunctions(junctions2);

	//	splittingRoad->initCalculation(numberOfTimesteps);
	//splittingRoad->loopTroughTimesteps();
	//	splittingRoad->dispResults();
	//	cout << endl << endl;
	//}




	////crossRoads
	//{

	//	OneWayRoadSSJ* crossRoads = new OneWayRoadSSJ(31, maxVelocity, dawdlePossibility, vehicleLength);


	//	vector<Source*> sources;
	//	sources.push_back(new SourceTerm(0, 1, crossRoads));
	//	sources.push_back(new SourceTerm(11, 1, crossRoads));

	//	crossRoads->setSources(sources);

	//	vector<Sink*> sinks3;
	//	sinks3.push_back(new SinkSimple(30, 0));
	//	sinks3.push_back(new SinkSimple(25, 0));
	//	crossRoads->setSinks(sinks3);

	//	vector<unsigned int> in3 = { 9,19 };
	//	vector<unsigned int> out3 = { 21,26 };
	//	vector<Junction*> junctions3 = { new JunctionSimple(in3, out3, crossRoads) };
	//	crossRoads->setJunctions(junctions3);

	//	crossRoads->initCalculation(numberOfTimesteps);
	//crossRoads->loopTroughTimesteps();
	//	crossRoads->dispResults();
	//	cout << endl << endl;
	//}



	//junctionTest(maxVelocity, vehicleLength);

	//junctionTestCrossRoads(maxVelocity, vehicleLength);

	std::cin.get();
}

