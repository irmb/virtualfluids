#ifndef AN_RESULT_TAYLOR_GREEN_VORTEX_H
#define AN_RESULT_TAYLOR_GREEN_VORTEX_H

#include "../AnalyticalResultProvider.h"
#include "VirtualFluids_GPU/LBM/LB.h"

class Results;
class TaylorGreenTestCondition;

class TaylorGreenAnalytical : public AnalyticalResultProvider
{
public:
	TaylorGreenAnalytical(std::vector< std::shared_ptr<Results> > simulationResults,std::shared_ptr<TaylorGreenTestCondition> testCondition);
	std::vector < std::shared_ptr<Results> > getAnalyticalResults();

private:
	void init(std::vector< std::shared_ptr<Results> > simulationResults);
	void calculate();

	real t;
	real x, z;
	real vx, vz;
	real press;
	real Amp;
	real L, L0, Lx, Lz;
	real rho0, vis;
	
};
#endif