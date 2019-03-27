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
	virtual void calculateTimeStep(TrafficMovement &road);

	virtual const std::vector<uint>& getInCellIndices() const;
	virtual const std::vector<uint>& getOutCellIndices() const;
	virtual const std::vector<bool>& getCarCanEnter() const;
	virtual const std::vector<int>& getCarsOnJunction()const;
	virtual const std::vector<uint>& getAlreadyMoved()const;
	virtual const std::vector<uint>& getOldSpeeds()const;
	virtual const std::vector<int>& getCarCanNotEnterThisOutCell()const;
	
	virtual void dispJunction(const uint index, const uint roadLength) const;
	virtual uint getNumCarsOnJunction() const; 

	virtual void checkOutCellIndices(const uint roadLength) const; 

private:
	uint getInCellsVectorIndex(uint cellIndex);

	void applyRules(int &carSpeed,int index, TrafficMovement &road);
	void brakeCar(uint outCellIndex, int &speed, int &remainingDistance, TrafficMovement &road);
	void moveCar(uint outCell, int carSpeed, int index, int remainingDistance, TrafficMovement &road);
	int chooseOutCell(int index);
	int generateRandomOutCellIndex(uint outCellsTempSize);

	void writeConcentrations(TrafficMovement &road);

};

