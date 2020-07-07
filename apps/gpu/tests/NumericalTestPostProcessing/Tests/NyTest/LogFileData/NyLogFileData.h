#ifndef NY_LOG_FILE_DATA_H
#define NY_LOG_FILE_DATA_H

#include <memory>
#include <string>
#include <vector>

class NyLogFileData
{
public:
	virtual std::vector<double> getBasicGridLengths() = 0;
	virtual int getStartTimeStepCalculation() = 0;
	virtual int getEndTimeStepCalculation() = 0;
	virtual std::string getDataToCalc() = 0;
	virtual std::vector<double> getNy() = 0;
	virtual std::vector<double> getNyDiff() = 0;
	virtual std::vector<std::vector<double> > getOrderOfAccuracy() = 0;
};
#endif