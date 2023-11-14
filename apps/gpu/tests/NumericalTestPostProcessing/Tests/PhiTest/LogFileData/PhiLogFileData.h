#ifndef PHI_LOG_FILE_DATA_H
#define PHI_LOG_FILE_DATA_H

#include <memory>
#include <string>
#include <vector>

class PhiLogFileData
{
public:
    virtual std::vector<double> getBasicGridLengths() = 0;
    virtual int getStartTimeStepCalculation() = 0;
    virtual int getEndTimeStepCalculation() = 0;
    virtual std::string getDataToCalc() = 0;
    virtual std::vector<double> getPhiDiff() = 0;
    virtual std::vector<std::vector<double> > getOrderOfAccuracy() = 0;
};
#endif