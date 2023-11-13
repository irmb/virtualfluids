#ifndef NY_TEST_PARAMETER_STRUCT_H
#define NY_TEST_PARAMETER_STRUCT_H

#include <memory>

#include "Utilities/Structs/BasicTestParameterStruct.h"

struct NyTestParameterStruct
{
    std::shared_ptr<BasicTestParameterStruct> basicTestParameter;

    double minOrderOfAccuracy;
    unsigned int startTimeStepCalculation;
    unsigned int endTimeStepCalculation;
};
#endif 