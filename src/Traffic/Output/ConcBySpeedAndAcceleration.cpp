#include "ConcBySpeedAndAcceleration.h"

#include <iostream>

ConcBySpeedAndAcceleration::ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed)
{
	std::cout << "using ConcBySpeedAndAcceleration::concentration for concentrations" << std::endl;
	concentration.resize(roadlength);
	this->maxSpeed = static_cast<real>(maxSpeed);
}


ConcBySpeedAndAcceleration::ConcBySpeedAndAcceleration(uint roadlength, real * concArrayStart, uint maxSpeed)
{
	if (concArrayStart == nullptr) {
		std::cout << "using ConcBySpeedAndAcceleration::concentration for concentrations" << std::endl;
		concentration.resize(roadlength);
		this->maxSpeed = static_cast<real>(maxSpeed);
	}
	else {
		std::cout << "using passed array for concentrations" << std::endl;
		useLBMConcArray = true;
		this->roadLength = roadlength;
		this->concArrayStart = concArrayStart;
		this->maxSpeed = static_cast<real>(maxSpeed);
	}

}

//
//void ConcBySpeedAndAcceleration::calculateConcForSingleCar(uint index, DrivingStates state, uint speed, uint acceleration)
//{	
//	if (useLBMConcArray) {
//		real *pos = concArrayStart + index;
//		*pos = static_cast<real>(speed) / maxSpeed;
//	}
//	else
//		concentration[index] = static_cast<real>(speed) / maxSpeed;	
//}


void ConcBySpeedAndAcceleration::calculateConcForSingleCar(uint index, uint oldSpeed, uint speed)
{
	putConcIntoArrayOrVector(index, chooseConc(oldSpeed, speed));
}


void ConcBySpeedAndAcceleration::calculateConcForJunctionCar(uint index, uint oldSpeed, uint speed)
{
	addConcToArrayOrVector(index, chooseConc(oldSpeed, speed));

}

real ConcBySpeedAndAcceleration::chooseConc(uint oldSpeed, uint speed)
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
