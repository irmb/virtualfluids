#include <iostream>
#include <vector>
#include <memory>

#include "Traffic/RoadMaker.h"
#include "Traffic/OneWayRoad.h"
#include "Traffic/OneWayRoadSSJ.h"
#include "Traffic/SourceTerm.h"
#include "Traffic/Junction.h"
#include "Traffic/JunctionSimple.h"
#include "Traffic/Sink.h"
#include "Traffic/SinkSimple.h"
#include "Traffic/ConcentrationByPosition.h"


using namespace std;

//static void junctionTest(unsigned int maxVelocity, unsigned int vehicleLength) {
//	{
//		float dawdlePossibility = 0.5;
//		vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, 2, -1, -1,-1,0,-1,-1,-1,-1,1,-1,-1,-1 };
//
//		OneWayRoadSSJ* junctionTest = new OneWayRoadSSJ(fiveCars, maxVelocity, dawdlePossibility, vehicleLength);
//
//		vector <unsigned int> in4 = { 9 };
//		vector<unsigned int> out4 = { 10 };
//		vector<Junction*> junctions4 = { new JunctionSimple(in4, out4, junctionTest) };
//		junctionTest->setJunctions(junctions4);
//
//		junctionTest->calculateResults(10);
//		junctionTest->dispResults();
//
//		cout << "Number of Cars" << junctionTest->getNumberOfCars() << endl;
//	}
//}
//
//static void junctionTestCrossRoads(unsigned int maxVelocity, unsigned int vehicleLength) {
//	{
//		float dawdlePossibility = 0.5;
//		vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, 2, -1, -1,-1,0,-1,-1,-1,-1,1,-1,-1,-1 };
//		OneWayRoadSSJ* junctionTestCrossRoads = new OneWayRoadSSJ(fiveCars, maxVelocity, dawdlePossibility, vehicleLength);
//
//		vector <unsigned int> in5 = { 6,9 };
//		vector<unsigned int> out5 = { 10, 14 };
//		vector<Junction*> junctions5 = {new JunctionSimple(in5, out5, junctionTestCrossRoads) };
//		junctionTestCrossRoads->setJunctions(junctions5);
//		junctionTestCrossRoads->setNeighbor(13, 7);
//
//		junctionTestCrossRoads->calculateResults(10);
//		junctionTestCrossRoads->dispResults();
//
//		cout << "Number of Cars" << junctionTestCrossRoads->getNumberOfCars() << endl;
//	}
//}


int main()
{
	int roadLength = 20;
	const unsigned int maxVelocity = 5;
	float vehicleDensity = 0.5;
	float dawdlePossibility = 0.5;//typical value: 0.2
	int numberOfTimesteps = 25;
	unsigned int vehicleLength = 3;
	vector<int> initialDistribution = { 1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,1,-1,-1,-1,-1 };
	vector<int> oneCar = { 1,-1,-1,-1,-1,-1 };
	vector<int> fiveCars = { 1, -1,-1, 1, -1, -1, -1, 2, -1, -1,-1,0,-1,-1,-1,-1,1,-1,-1,-1 };

	//OneWay random 
	{
		cout << "OneWay random" << endl;


		auto r = make_shared<RoadMaker>(roadLength, maxVelocity, vehicleLength, vehicleDensity);
		auto road = make_shared<OneWayRoad>(r ,dawdlePossibility);

		//auto writer = make_unique<ConcentrationByPosition>(ConcentrationByPosition(road->getRoadLength()));
		//road->setConcentrationOutwriter(move(writer));

		road->calculateResults(numberOfTimesteps);
		road->dispResults();

		cout << endl << endl;

	}


	////OneWay initial distribution
	//{
	//	cout << "OneWay initial distribution" << endl;
	//	unique_ptr<RoadNetworkData> r = make_unique<RoadMaker>(initialDistribution, maxVelocity, vehicleLength);
	//	shared_ptr<OneWayRoad> road = make_shared<OneWayRoad>(move(r), dawdlePossibility);
	//	road->calculateResults(numberOfTimesteps);
	//	road->dispResults();	
	//	cout << endl << endl;
	//}


	//////OneWay slowToStart
	//{
	//	dawdlePossibility = 0;
	//	cout << "OneWay slowToStart" << endl;
	//	vector<int> initialDist = { -1,1,5,-1,-1,-1,-1,-1,-1,1,-1,-1,-1,0,-1,-1,-1,-1,-1 };

	//	unique_ptr<RoadNetworkData> r = make_unique<RoadMaker>(initialDist, maxVelocity, 1);
	//	shared_ptr<OneWayRoad> road = make_shared<OneWayRoad>(move(r), dawdlePossibility);
	//	road->setSlowToStart(static_cast<float>(0.9999));
	//	road->calculateResults(numberOfTimesteps);
	//	road->dispResults();
	//	cout << endl << endl;
	//}



	////sources and sinks
	//{
	//cout << "sources and sinks" << endl;

	//vector< unique_ptr<Sink> > sinks;
	//sinks.push_back( make_unique <SinkSimple>(5, static_cast<float>(0.5)));

	//unique_ptr<RoadMaker> r = make_unique<RoadMaker>(oneCar, maxVelocity, vehicleLength);

	//vector< unique_ptr<Source> > sources;
	//sources.push_back( make_unique <SourceTerm>(0, static_cast<float>(1.0), r->getMaxVelocity()) );

	//r->setSources(sources);
	//r->setSinks(sinks);

	//shared_ptr<OneWayRoadSSJ> roadSource = make_shared<OneWayRoadSSJ>(move(r), static_cast<float>(0.0));

	//roadSource->calculateResults(numberOfTimesteps);
	//roadSource->dispResults();
	//cout << endl << endl;
	//}



	////mergingRoad
	//{
	//	unique_ptr<RoadMaker> r = make_unique<RoadMaker>(25, maxVelocity, vehicleLength);
	//	
	//	vector< unique_ptr<Source> > sources;
	//	sources.push_back(make_unique <SourceTerm>(0, static_cast<float>(0.5), r->getMaxVelocity()));
	//	sources.push_back(make_unique <SourceTerm>(10, static_cast<float>(0.5), r->getMaxVelocity()));
	//	r->setSources(sources);

	//	unique_ptr<Sink> s = make_unique <SinkSimple>(SinkSimple(24, static_cast<float>(0.0)));
	//	r->addSink(s);

	//	vector< unique_ptr<Junction> > junctions;
	//	vector<unsigned int> in = { 9,20 };
	//	vector<unsigned int> out = { 21 };
	//	junctions.push_back(make_unique <JunctionSimple>(JunctionSimple(in, out)));

	//	r->setJunctions(junctions);

	//	shared_ptr<OneWayRoadSSJ> mergingRoad = make_shared<OneWayRoadSSJ>(move(r), dawdlePossibility);


	//	mergingRoad->calculateResults(numberOfTimesteps);
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

	//	splittingRoad->calculateResults(numberOfTimesteps);
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

	//	crossRoads->calculateResults(numberOfTimesteps);
	//	crossRoads->dispResults();
	//	cout << endl << endl;
	//}



	//junctionTest(maxVelocity, vehicleLength);

	//junctionTestCrossRoads(maxVelocity, vehicleLength);

	std::cin.get();
}

