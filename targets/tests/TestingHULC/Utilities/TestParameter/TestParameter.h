#ifndef TESTPARAMETERINT_H
#define TESTPARAMETERINT_H

#include <memory>
#include <string>

class InitialCondition;

class TestParameter
{
public:
	virtual std::shared_ptr<InitialCondition> getInitialCondition() = 0;
	virtual double getViscosity() = 0;
	virtual std::string getGridPath() = 0;
	virtual std::string getFilePath() = 0;
	virtual unsigned int getNumberOfGridLevels() = 0;
	virtual unsigned int getEndTime() = 0;
	virtual unsigned int getTimeStepLength() = 0;
	virtual unsigned int getLx() = 0;
	virtual unsigned int getLz() = 0;
	virtual unsigned int getYSliceForCalculation() = 0;
	virtual unsigned int getStartTimeCalculation() = 0;
	virtual bool getWriteFiles() = 0;
	virtual unsigned int getStartTimeDataWriter() = 0;

private:

};

#endif // !TESTPARAMETER_H
