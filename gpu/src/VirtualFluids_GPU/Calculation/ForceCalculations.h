#ifndef FORCE_CALCULATIONS_H
#define FORCE_CALCULATIONS_H

#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"

#include <iostream>
#include <stdio.h>

class CudaMemoryManager;

class ForceCalculations
{
public:
	ForceCalculations(Parameter* para);
	~ForceCalculations();
	void calcPIDControllerForForce(Parameter* para, CudaMemoryManager* cudaManager);
	void printForcing(Parameter* para);

private:
	double vx1Targed; //!< target velocity.
	double Kpcrit; //Kp critical
	double Tcrit;  //the oscillation period 
	double Tn;
	double Tv;
	double e;
	double Ta;
	double Kp;
	double Ki;
	double Kd;
	double y;
	double esum;
	double eold;
	bool isPID;
};


//extern "C" void calcVeloForce(Parameter* para);
//extern "C" void allocVeloForForcing(Parameter* para);
//extern "C" void printForcing(Parameter* para);

#endif /* FORCE_CALCULATIONS_H */
