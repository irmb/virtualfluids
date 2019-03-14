#pragma once
#include <random> 
#include <vector>

#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"

#include "Junction.h"

class TrafficMovement;

class VF_PUBLIC JunctionRandom :
	public Junction
{

private:
	JunctionData data;

public:
	JunctionRandom(const std::vector<uint> &inCellIndices, const std::vector<uint> &outCellIndices);
	~JunctionRandom() {};

	virtual void setCellIndexForNoUTurn(std::vector<int> carCanNotEnterThisOutCell);

	virtual bool acceptsCar(uint cellIndex); //determines if a car can enter the junction
	virtual void registerCar(uint cellIndex, uint numberOfCellsAlreadyMoved,  uint speed, uint oldSpeed); //registers all cars entering the junction
	virtual void calculateTimeStep(TrafficMovement& road);
	virtual void updateJunction();

	virtual const std::vector<uint>& getInCellIndices() const;
	
	virtual void dispJunction(const uint index, uint roadLength) const;
	virtual const uint getNumCarsOnJunction() const; 

	virtual void checkOutCellIndices(uint roadLength); 

private:
	uint getInCellsVectorIndex(uint cellIndex);

	void applyRules(int &carSpeed,const int &index, TrafficMovement& road);
	void breakCar(uint outCellIndex, int &speed, int &remainingDistance, const int & index, TrafficMovement& road);
	void moveCar(uint outCell, int & carSpeed, const int & index, int remainingDistance, TrafficMovement& road);
	int chooseOutCell(const int & index);
	int generateRandomOutCellIndex(uint outCellsTempSize);

	void writeConcentrations(TrafficMovement& road);

};

