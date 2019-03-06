 #pragma once
#include <random>
#include "RoadNetworkData.h"
#include "Utilities/RandomHelper.h"
#include <VirtualFluidsDefinitions.h>

struct VF_PUBLIC RoadMaker :
	public RoadNetworkData
{
public:
	RoadMaker(const uint roadLength, const uint maxVelocity, uint vehicleLength, const float vehicleDensity); //random vehicle Distribution
	RoadMaker(const std::vector<int> vehicleDistribution, const uint maxVelocity, uint vehicleLength); //given vehicle distribution
	RoadMaker(const uint roadLength, const uint maxVelocity, uint vehicleLength);//empty road

	~RoadMaker();

	void setJunctions( std::vector<std::unique_ptr<Junction> > & junctions); //max 999 junctions
	void addJunction(std::unique_ptr<Junction> & junction);
	void setSinks(std::vector< std::unique_ptr<Sink> > & sinks); //max 999 sinks
	void addSink(std::unique_ptr<Sink> & sink);
	void setSources(std::vector< std::unique_ptr<Source> > & sources);
	void addSource(std::unique_ptr<Source> & source);

	void setNeighbor(uint index, uint neighbor); // don't use it for setting sinks or junctions!

	uint getMaxVelocity();

private:
	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<float> distFloat{ 0.0, 1.0 };
	std::uniform_int_distribution<uint> distInt{ 0, maxVelocity };

private:
	void initNext();
	void initNeighbors();
	void initCurrentAsEmpty();
	void initCurrentWithLongVehicles();
	void initVehicleDensity(const float vehicleDensity);
	void initRandomCars(const float vehicleDensity);
	void initVehicleLength(const uint vehicleLength);
	int randomSpeed();

	void setJunctionAsNeighbor(std::unique_ptr<Junction> & junction);
	void setSinkAsNeighbor(std::unique_ptr<Sink> & sink);

};

