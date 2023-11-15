#ifndef INITIAL_CONDITION_SHEAR_WAVE_H
#define    INITIAL_CONDITION_SHEAR_WAVE_H

#include "Utilities/InitialCondition/InitialConditionImp.h"

#include <memory>

struct ShearWaveParameterStruct;
struct GridInformationStruct;

class InitialConditionShearWave :public InitialConditionImp
{
public:
    static std::shared_ptr<InitialConditionShearWave> getNewInstance(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct);

    real getInitVX(int i, int level);
    real getInitVY(int i, int level);
    real getInitVZ(int i, int level);
    real getInitROH(int i, int level);
    real getInitPRESS(int i, int level);

private:
    InitialConditionShearWave(std::shared_ptr<ShearWaveParameterStruct> simParaStruct, std::shared_ptr<GridInformationStruct> gridInfoStruct);
    InitialConditionShearWave() {};
    
    real rho;
    real l0;
    real lx, lz;
    real u0, v0;
};
#endif 