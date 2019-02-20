#ifndef INITIAL_CONDITION_TAYLORGREENVORTEX_UX_H
#define INITIAL_CONDITION_TAYLORGREENVORTEX_UX_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

struct TaylorGreenVortexUxParameterStruct;
struct GridInformationStruct;

class InitialConditionTaylorGreenUx : public InitialConditionImp
{
public:
	static std::shared_ptr<InitialConditionTaylorGreenUx> getNewInstance(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct);

	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionTaylorGreenUx(std::shared_ptr<TaylorGreenVortexUxParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct>  gridInfoStruct);
	InitialConditionTaylorGreenUx() {};

	real amp;
	real rho;
	real l0;
	real lx, lz;
	real ux;
};

#endif