#pragma once
#include <random>
#include "RoadNetworkData.h"
#include "Utilities/RandomHelper.h"
#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC RoadMaker :
	public RoadNetworkData
{
public:
	RoadMaker(const unsigned int roadLength, const unsigned int maxVelocity, unsigned int vehicleLength, const float vehicleDensity); //random vehicle Distribution
	RoadMaker(const std::vector<int> vehicleDistribution, const unsigned int maxVelocity, unsigned int vehicleLength); //given vehicle distribution
	RoadMaker(const unsigned int roadLength, const unsigned int maxVelocity, unsigned int vehicleLength);//empty road

	~RoadMaker();

	void setJunctions( std::vector<std::unique_ptr<Junction> > & junctions); //max 999 junctions
	void addJunction(std::unique_ptr<Junction> & junction);
	void setSinks(std::vector< std::unique_ptr<Sink> > & sinks); //max 999 sinks
	void addSink(std::unique_ptr<Sink> & sink);
	void setSources(std::vector< std::unique_ptr<Source> > & sources);
	void addSource(std::unique_ptr<Source> & source);

	void setNeighbor(unsigned int index, unsigned int neighbor); // don't use it for setting sinks or junctions!

	unsigned int getMaxVelocity();

private:
	mt19937 engine = RandomHelper::make_engine();
	uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
	uniform_int_distribution<unsigned int> distInt{ 0, maxVelocity };

private:
	void initNext();
	void initNeighbors();
	void initCurrentAsEmpty();
	void initCurrentWithLongVehicles();
	void initVehicleDensity(const float vehicleDensity);
	void initRandomCars(const float vehicleDensity);
	void initVehicleLength(const unsigned int vehicleLength);
	int randomSpeed();

	void setJunctionAsNeighbor(std::unique_ptr<Junction> & junction);
	void setSinkAsNeighbor(std::unique_ptr<Sink> & sink);

};

