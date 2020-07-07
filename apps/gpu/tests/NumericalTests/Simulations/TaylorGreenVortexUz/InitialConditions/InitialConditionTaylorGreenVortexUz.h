#ifndef INITIAL_CONDITION_TAYLORGREENVORTEX_UZ_H
#define INITIAL_CONDITION_TAYLORGREENVORTEX_UZ_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

struct TaylorGreenVortexUzParameterStruct;
struct GridInformationStruct;

class InitialConditionTaylorGreenUz :public InitialConditionImp
{
public:
	static std::shared_ptr<InitialConditionTaylorGreenUz> getNewInstance(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct);

	real getInitVX(int i, int level);
	real getInitVY(int i, int level);
	real getInitVZ(int i, int level);
	real getInitROH(int i, int level);
	real getInitPRESS(int i, int level);

private:
	InitialConditionTaylorGreenUz(std::shared_ptr<TaylorGreenVortexUzParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct);
	InitialConditionTaylorGreenUz() {};

	real amp;
	real rho;
	real l0;
	real lx, lz;
	real uz;
};

#endif