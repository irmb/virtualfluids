#ifndef PHI_TEST_PARAMETER_STRUCT_H
#define PHI_TEST_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities/Structs/BasicTestParameterStruct.h"

struct PhiTestParameterStruct
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

    double minOrderOfAccuracy;
    unsigned int startTimeStepCalculation;
    unsigned int endTimeStepCalculation;
};
#endif 