#pragma once
#include <VirtualFluidsDefinitions.h>

#include <vector>
#include <iostream>
#include <random> 
#include <algorithm>

#include "Utilities/invalidInput_error.h"
#include "Utilities/VectorHelper.h"
#include "Junction.h"
#include "OneWayRoadSSJ.h"

using namespace std;

class VF_PUBLIC JunctionSimple :
	public Junction
{

private:
	JunctionData data;

public:
	JunctionSimple(const vector<unsigned int> &inCellIndices, const vector<unsigned int> &outCellIndices);
	~JunctionSimple() {};

	virtual void setCellIndexForNoUTurn(vector<int> carCanNotEnterThisOutCell);

	virtual bool acceptsCar(unsigned int cellIndex); //determines if a car can enter the junction
	virtual void registerCar(unsigned int cellIndex, unsigned int numberOfCellsAlreadyMoved,  unsigned int speed); //registers all cars entering the junction
	virtual void calculateTimeStep(OneWayRoadSSJ& road);
	virtual void updateJunction();

	virtual const vector<unsigned int>& getInCellIndices() const;
	
	virtual void dispJunction(const unsigned int index, unsigned int roadLength) const;
	virtual const unsigned int getNumCarsOnJunction() const; 

	virtual void checkOutCellIndices(unsigned int roadLength); 

private:
	unsigned int getInCellsVectorIndex(unsigned int cellIndex);

	void applyRules(int &carSpeed,const int &index, OneWayRoadSSJ& road);
	void breakCar(unsigned int outCellIndex, int &speed, unsigned int &remainingDistance, OneWayRoadSSJ& road);
	void moveCar(unsigned int outCell, int & carSpeed, const int & index, OneWayRoadSSJ& road);
	int chooseOutCell(const int & index);


private:
	//variables for temporaray calculations
	unsigned int remainingDistance;
	int outCell;
	unsigned int random;
	unsigned int gap;
	int index;
};

