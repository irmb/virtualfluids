#ifndef INITIAL_CONDITION_IMP_H
#define INITIAL_CONDITION_IMP_H

#include "InitialCondition.h"

#include "VirtualFluids_GPU/LBM/LB.h"

#include <vector>
#include <memory>

class Parameter;

class InitialConditionImp : public InitialCondition
{
public:
	void setParameter(std::shared_ptr<Parameter> para);
	void init(const int level);
	virtual real getInitVX(int i, int level) = 0;
	virtual real getInitVY(int i, int level) = 0;
	virtual real getInitVZ(int i, int level) = 0;
	virtual real getInitROH(int i, int level) = 0;
	virtual real getInitPRESS(int i, int level) = 0;

protected:
	InitialConditionImp() {};
	real getXCoord(int i, int level);
	real getYCoord(int i, int level);
	real getZCoord(int i, int level);

	std::shared_ptr<Parameter> para;
	real XCoordstopnode, YCoordstopnode, ZCoordstopnode;

};
#endif