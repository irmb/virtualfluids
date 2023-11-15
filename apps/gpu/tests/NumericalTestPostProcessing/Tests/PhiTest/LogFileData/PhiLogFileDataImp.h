#ifndef PHI_LOG_FILE_DATA_IMP_H
#define PHI_LOG_FILE_DATA_IMP_H

#include "PhiLogFileData.h"

#include <memory>

class PhiLogFileDataImp : public PhiLogFileData
{
public:
    static std::shared_ptr<PhiLogFileDataImp> getNewInstance();

    std::vector<double> getBasicGridLengths();
    int getStartTimeStepCalculation();
    int getEndTimeStepCalculation();
    std::string getDataToCalc();
    std::vector<double> getPhiDiff();
    std::vector<std::vector<double> > getOrderOfAccuracy();

    void setBasicGridLengths(std::vector<double> basicGridLengths);
    void setStartTimeStepCalculation(int startTimeStepCalculation);
    void setEndTimeStepCalculation(int endTimeStepCalculation);
    void setDataToCalc(std::string dataToCalcPhiAndNu);
    void setPhiDiff(std::vector<double> phiDiff);
    void setOrderOfAccuracy(std::vector<std::vector<double> > orderOfAccuracy);

private:
    PhiLogFileDataImp();

    std::vector<double> basicGridLengths;
    int startTimeStepCalculation;
    int endTimeStepCalculation;
    std::string dataToCalc;
    std::vector<double> phiDiff;
    std::vector<std::vector<double> > orderOfAccuracy;
};
#endif