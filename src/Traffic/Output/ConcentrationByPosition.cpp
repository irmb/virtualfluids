#include "ConcentrationByPosition.h"



ConcentrationByPosition::ConcentrationByPosition(uint roadlength, uint maxSpeed, uint maxAcceleration)
{
	concentration.resize(roadlength);
}


//void ConcentrationByPosition::calculateConcFromCarDistribution(const std::vector<int>& currentCarDistribution)
//{
//	for (uint i = 0; i < currentCarDistribution.size(); i++) {
//		if (currentCarDistribution[i] >= 0) 
//			concentration[i] = 1.0;
//		else
//			concentration[i] = 0.0;
//	}
//
//	//dispConcentration();
//}


void ConcentrationByPosition::calculateConcForSingleCar(uint index, uint speed, uint acceleration)
{
	concentration[index] = 1.0;
}



