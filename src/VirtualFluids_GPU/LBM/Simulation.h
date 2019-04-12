#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <memory>
#include "Core/PointerDefinitions.h"

#include "VirtualFluidsDefinitions.h"

#include "Output/LogWriter.hpp"
#include "GPU/KineticEnergyAnalyzer.h"
#include "GPU/EnstrophyAnalyzer.h"
#include "Utilities/Buffer2D.hpp"
#include "LBM/LB.h"

class Communicator;
class Parameter;
class GridProvider;
class PorousMedia;
class RestartObject;
class RestartPostprocessor;
class ForceCalculations;
class DataWriter;
class Kernel;
class TrafficMovementFactory;


class VF_PUBLIC Simulation
{
public:
	Simulation();
	~Simulation();
	void run();
	void init(SPtr<Parameter> para, SPtr<GridProvider> gridProvider, std::shared_ptr<DataWriter> dataWriter);
	void free();
	void bulk();
	void porousMedia();
	void definePMarea(PorousMedia* pm);

protected:
	Buffer2D <real> sbuf_t; 
	Buffer2D <real> rbuf_t;
	Buffer2D <real> sbuf_b;
	Buffer2D <real> rbuf_b;

	Buffer2D <int> geo_sbuf_t; 
	Buffer2D <int> geo_rbuf_t;
	Buffer2D <int> geo_sbuf_b;
	Buffer2D <int> geo_rbuf_b;


	LogWriter output;

    Communicator* comm;
    SPtr<Parameter> para;
    SPtr<GridProvider> gridProvider;
    SPtr<DataWriter> dataWriter;
	std::vector < SPtr< Kernel>> kernels;

	//Restart object
	RestartObject* restObj;
	RestartPostprocessor* rest;

	//Forcing Calculation
	ForceCalculations* forceCalculator;

	//Porous Media
	PorousMedia* pm[3];
	//PorousMedia* pm0;
	//PorousMedia* pm1;
	//PorousMedia* pm2;

	//KQ - Schlaff
	unsigned int            kNQ, kSQ, kEQ, kWQ;
	QforBoundaryConditions  QnH, QnD;
	QforBoundaryConditions  QsH, QsD;
	QforBoundaryConditions  QeH, QeD;
	QforBoundaryConditions  QwH, QwD;
	real *VxNH,          *VyNH,       *VzNH,       *deltaVNH;
	real *VxND,          *VyND,       *VzND,       *deltaVND;
	real *VxSH,          *VySH,       *VzSH,       *deltaVSH;
	real *VxSD,          *VySD,       *VzSD,       *deltaVSD;
	real *VxEH,          *VyEH,       *VzEH,       *deltaVEH;
	real *VxED,          *VyED,       *VzED,       *deltaVED;
	real *VxWH,          *VyWH,       *VzWH,       *deltaVWH;
	real *VxWD,          *VyWD,       *VzWD,       *deltaVWD;


	////////////////////////////////////////////////////////////////////////////
	KineticEnergyAnalyzer *kinEnergyWriter;
	////////////////////////////////////////////////////////////////////////////
	EnstrophyAnalyzer *enstrophyWriter;
	////////////////////////////////////////////////////////////////////////////
	TrafficMovementFactory *trafficFactory;
	bool useTrafficGPU;

 };
#endif
