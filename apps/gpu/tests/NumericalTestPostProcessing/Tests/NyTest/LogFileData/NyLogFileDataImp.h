#ifndef NY_LOG_FILE_DATA_IMP_H
#define NY_LOG_FILE_DATA_IMP_H

#include "NyLogFileData.h"

#include <memory>

class NyLogFileDataImp : public NyLogFileData
{
public:
    static std::shared_ptr<NyLogFileDataImp> getNewInstance();

    std::vector<double> getBasicGridLengths();
    int getStartTimeStepCalculation();
    int getEndTimeStepCalculation();
    std::string getDataToCalc();
    std::vector<double> getNy();
    std::vector<double> getNyDiff();
    std::vector<std::vector<double> > getOrderOfAccuracy();

    void setBasicGridLengths(std::vector<double> basicGridLengths);
    void setStartTimeStepCalculation(int startTimeStepCalculation);
    void setEndTimeStepCalculation(int endTimeStepCalculation);
    void setDataToCalc(std::string dataToCalcPhiAndNu);
    void setNy(std::vector<double> ny);
    void setNyDiff(std::vector<double> nyDiff);
    void setOrderOfAccuracy(std::vector<std::vector<double> > orderOfAccuracy);


private:
    NyLogFileDataImp();

    std::vector<double> basicGridLengths;
    int startTimeStepCalculation;
    int endTimeStepCalculation;
    std::string dataToCalc;
    std::vector<double> ny;
    std::vector<double> nyDiff;
    std::vector<std::vector<double> > orderOfAccuracy;

};
#endif