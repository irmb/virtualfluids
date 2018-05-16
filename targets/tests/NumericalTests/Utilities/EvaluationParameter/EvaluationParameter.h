#ifndef EVALUATIONPARAMETER_H
#define EVALUATIONPARAMETER_H

#include <memory>
#include <string>
#include <time.h>

class EvaluationParameter
{
public:
	static std::shared_ptr<EvaluationParameter> getNewInstance(std::string testNames, int numberOfEqualTests, int lx, std::string dataToCalculate, std::string logFilePath, double minOrderOfAccuracy, bool writeFiles, double viscosity);
	std::string getTestName();
	std::string getDataToCalculate();
	int getLx();
	int getNumberOfEqualTests();
	void setStartTime();
	void setEndTime();
	double getTestTime();
	std::string getLogFilePath();
	double getMinOrderOfAccuracy();
	bool getWriteFiles();
	double getViscosity();

protected:
	EvaluationParameter() {};
	EvaluationParameter(std::string testNames, int numberOfEqualTests, int lx, std::string dataToCalculate, std::string logFilePath, double minOrderOfAccuracy, bool writeFiles, double viscosity);

private:
	int lx;
	std::string testName;
	int numberOfEqualTests;
	std::string dataToCalculate;
	time_t startTime, endTime;
	std::string logFilePath;
	double minOrderOfAccuracy;
	bool writeFiles;
	double viscosity;
};
#endif 
