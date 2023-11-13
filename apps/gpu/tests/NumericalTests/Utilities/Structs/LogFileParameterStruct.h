#ifndef LOG_FILE_PARAMETER_STRUCT_H
#define LOG_FILE_PARAMETER_STRUCT_H

#include <vector>

struct LogFileParameterStruct
{
    std::vector<int> devices;
    int numberOfTimeSteps;
    bool writeAnalyticalToVTK;

};
#endif 