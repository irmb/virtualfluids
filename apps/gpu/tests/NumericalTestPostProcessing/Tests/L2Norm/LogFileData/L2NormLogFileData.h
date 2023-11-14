#ifndef L2_NORM_LOG_FILE_DATA_H
#define L2_NORM_LOG_FILE_DATA_H

#include <string>
#include <vector>

class L2NormLogFileData
{
public:
    virtual std::vector<double> getBasicGridLengths() = 0;
    virtual std::string getDataToCalc() = 0;
    virtual std::string getNormalizeData() = 0;
    virtual int getBasicTimeStep() = 0;
    virtual int getDivergentTimeStep() = 0;
    virtual std::vector<double> getL2NormForBasicTimeStep() = 0;
    virtual std::vector<double> getL2NormForDivergentTimeStep() = 0;
    virtual std::vector<double> getL2NormDiff() = 0;
};
#endif