#include "ConcBySpeedAndAcceleration.h"

ConcBySpeedAndAcceleration::ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed)
{
	concentration.resize(roadlength);
	this->maxAcceleration = static_cast<float>(maxAcceleration);
	this->maxSpeed = static_cast<float>(maxSpeed);
}


ConcBySpeedAndAcceleration::ConcBySpeedAndAcceleration(uint roadlength, float * concArrayStart, uint concArraySize, uint maxSpeed)
{
	useLBMConcArray = true;
	this->concArraySize = concArraySize;
	this->concArrayStart = concArrayStart;
	this->maxAcceleration = static_cast<float>(maxAcceleration);
	this->maxSpeed = static_cast<float>(maxSpeed);
}

//
//void ConcBySpeedAndAcceleration::calculateConcForSingleCar(uint index, DrivingStates state, uint speed, uint acceleration)
//{	
//	if (useLBMConcArray) {
//		float *pos = concArrayStart + index;
//		*pos = static_cast<float>(speed) / maxSpeed;
//	}
//	else
//		concentration[index] = static_cast<float>(speed) / maxSpeed;	
//}


void ConcBySpeedAndAcceleration::calculateConcForSingleCar(uint index, uint oldSpeed, uint speed)
{
	putConcIntoArrayOrVector(index, chooseConc(oldSpeed, speed));
}


void ConcBySpeedAndAcceleration::calculateConcForJunctionCar(uint index, uint oldSpeed, uint speed)
{
	addConcToArrayOrVector(index, chooseConc(oldSpeed, speed));

}

float ConcBySpeedAndAcceleration::chooseConc(uint oldSpeed, uint speed)
{
	if (oldSpeed == 0 && speed > 0) //Start
		return 1.0f;
	else if (oldSpeed == 0 && speed == 0) //Idle
		return 0.25f;
	else if (speed == oldSpeed) //Drive
		return 0.5f;
	else if (speed > oldSpeed) //Accelerate
		return 0.75f;
	else if (speed < oldSpeed) //Brake
		return 0.45f;
	else 
		std::cerr << "couldn't choose driving state in ConcentrationBySpeedAndAcceleration::chooseConc" << std::endl;
	return -1.0f;
}
