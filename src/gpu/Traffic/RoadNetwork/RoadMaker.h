 #pragma once
#include <random>

#include "RoadNetworkData.h"

#include "Utilities/RandomHelper.h"


struct VIRTUALFLUIDS_GPU_EXPORT RoadMaker :
	public RoadNetworkData
{
public:
	RoadMaker(const uint roadLength, const uint maxVelocity, uint vehicleLength, const real vehicleDensity); //random vehicle Distribution
	RoadMaker(const std::vector<int> vehicleDistribution, const uint maxVelocity, uint vehicleLength); //given vehicle distribution
	RoadMaker(const uint roadLength, const uint maxVelocity, uint vehicleLength);//empty road

	~RoadMaker();

	void setJunctions( std::vector<std::shared_ptr<Junction> > & junctions); //max 999 junctions
	void addJunction(std::shared_ptr<Junction> & junction);
	void setSinks(std::vector< std::shared_ptr<Sink> > & sinks); //max 999 sinks
	void addSink(std::shared_ptr<Sink> & sink);
	void setSources(std::vector< std::shared_ptr<Source> > & sources);
	void addSource(std::shared_ptr<Source> & source);

	void setNeighbor(uint index, uint neighbor); // don't use it for setting sinks or junctions!
	void setNeighborForCurve(uint index, uint neighbor);

	uint getMaxVelocity();

private:
	std::mt19937 engine = RandomHelper::make_engine();
	std::uniform_real_distribution<real> distFloat{ 0.0, 1.0 };
	std::uniform_int_distribution<uint> distInt{ 0, maxVelocity };

private:
	void initNext();
	void initNeighbors();
	void initCurrentAsEmpty();
	void initCurrentWithLongVehicles();
	void initVehicleDensity(const real vehicleDensity);
	void initRandomCars(const real vehicleDensity);
	void initVehicleLength(const uint vehicleLength);
	int randomSpeed();

	void setJunctionAsNeighbor(std::shared_ptr<Junction> & junction);
	void setSinkAsNeighbor(std::shared_ptr<Sink> & sink);

};

