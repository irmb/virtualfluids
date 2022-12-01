
#ifndef TESTUTILITIESGPU_H
#define TESTUTILITIESGPU_H

#include "../src/gpu/VirtualFluids_GPU/Parameter/Parameter.h"

namespace testingVF
{

inline SPtr<Parameter> createParameterForLevel(uint level)
{
    SPtr<Parameter> para = std::make_shared<Parameter>();
    para->setMaxLevel(level + 1); // setMaxLevel resizes parH and parD
    para->parH[level] = std::make_shared<LBMSimulationParameter>();
    para->parD[level] = std::make_shared<LBMSimulationParameter>();

    return para;
}

} // namespace testingVF

#endif