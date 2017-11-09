#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include "Output/LogWriter.hpp"
#include "Communication/Communicator.h"
#include "Utilities/Buffer2D.hpp"
//#include "Input/ConfigFile.h"
#include "Calculation/ForceCalculations.h"
#include "Parameter/Parameter.h"
#include "Restart/RestartPostprocessor.h"
#include "Restart/RestartObject.h"
#include "Utilities/StringUtil.hpp"
#include "LBM/LB.h"
//#include "LBM/D3Q27.h"


class Simulation
{
public:
	Simulation(void);
	~Simulation(void);
	void run();
	void init(std::string &cstr);
	void bulk();
protected:

	Buffer2D <doubflo> sbuf_t; 
	Buffer2D <doubflo> rbuf_t;
	Buffer2D <doubflo> sbuf_b;
	Buffer2D <doubflo> rbuf_b;

	Buffer2D <int> geo_sbuf_t; 
	Buffer2D <int> geo_rbuf_t;
	Buffer2D <int> geo_sbuf_b;
	Buffer2D <int> geo_rbuf_b;

	//for MPI
	Communicator* comm;
	LogWriter output;

	//Parameter
	Parameter* para;

	//Restart object
	RestartObject* restObj;
	RestartPostprocessor* rest;

	//Forcing Calculation
	ForceCalculations* forceCalculator;

	//KQ - Schlaff
	unsigned int            kNQ, kSQ, kEQ, kWQ;
	QforBoundaryConditions  QnH, QnD;
	QforBoundaryConditions  QsH, QsD;
	QforBoundaryConditions  QeH, QeD;
	QforBoundaryConditions  QwH, QwD;
	doubflo *VxNH,          *VyNH,       *VzNH,       *deltaVNH;
	doubflo *VxND,          *VyND,       *VzND,       *deltaVND;
	doubflo *VxSH,          *VySH,       *VzSH,       *deltaVSH;
	doubflo *VxSD,          *VySD,       *VzSD,       *deltaVSD;
	doubflo *VxEH,          *VyEH,       *VzEH,       *deltaVEH;
	doubflo *VxED,          *VyED,       *VzED,       *deltaVED;
	doubflo *VxWH,          *VyWH,       *VzWH,       *deltaVWH;
	doubflo *VxWD,          *VyWD,       *VzWD,       *deltaVWD;
 };
#endif
