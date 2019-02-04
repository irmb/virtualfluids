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

	real Amp;
	real rho;
	real L0;
	real Lx, Lz;
	real ux;
};

#endif