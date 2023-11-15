#ifndef L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_H
#define L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_H

#include <string>
#include <vector>

class L2NormBetweenKernelsLogFileData
{
public:
    virtual std::vector<double> getBasicGridLengths() = 0;
    virtual std::string getBasicKernel() = 0;
    virtual std::string getDivergentKernel() = 0;
    virtual std::string getDataToCalculate() = 0;
    virtual int getTimeStep() = 0;
    virtual std::vector<double> getL2NormForBasicKernel() = 0;
    virtual std::vector<double> getL2NormForDivergentKernel() = 0;
    virtual std::vector<double> getL2NormBetweenKernels() = 0;
    virtual std::string getNormalizeData() = 0;
};
#endif