#include "ConcBySpeedAndAcceleration.h"

ConcBySpeedAndAcceleration::ConcBySpeedAndAcceleration(uint roadlength, uint maxSpeed, uint maxAcceleration)
{
	concentration.resize(roadlength);
	this->maxAcceleration = static_cast<float>(maxAcceleration);
	this->maxSpeed = static_cast<float>(maxSpeed);
}


void ConcBySpeedAndAcceleration::calculateConcForSingleCar(uint index, uint speed, uint acceleration)
{
	concentration[index] = static_cast<float>(speed) / maxSpeed;
}