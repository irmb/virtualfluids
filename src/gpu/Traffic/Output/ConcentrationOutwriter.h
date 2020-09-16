#pragma once

#include <vector>


#include "Core/DataTypes.h"

#include "Traffic_export.h"

class TRAFFIC_EXPORT ConcentrationOutwriter
{
public:
	virtual void resetConcentrations();
	virtual void calculateConcForSingleCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;
	virtual void calculateConcForJunctionCar(uint index, uint oldSpeed = 0, uint speed = 0) = 0;
	virtual void calculateConcForAllCars(const std::vector<int> oldSpeeds, const std::vector<int> newSpeeds)=0;
	void dispCurrentConcentrations();

protected:
	void putConcIntoArrayOrVector(uint index, real conc);
	void addConcToArrayOrVector(uint index, real conc);

protected:
	std::vector<real> concentration;
	bool useLBMConcArray = false;
	real* concArrayStart;
	uint roadLength;

private:
	void dispSingleConcentration(real conc);
	
};

