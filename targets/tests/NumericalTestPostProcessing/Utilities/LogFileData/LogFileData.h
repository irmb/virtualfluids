#ifndef LOGFILE_DATA_H
#define LOGFILE_DATA_H

#include <memory>
#include <string>
#include <vector>

class L2NormLogFileData;
class L2NormBetweenKernelsLogFileData;
class NyLogFileData;
class PhiLogFileData;
class ShearWaveLogFileData;
class TaylorGreenVortexUxLogFileData;
class TaylorGreenVortexUzLogFileData;

class LogFileData
{
public:
	virtual std::vector<int> getAnalyticalVTKWritingTime() = 0;
	virtual std::vector<double> getBasicGridLengths() = 0;
	virtual int getBasisTimeStepLength() = 0;
	virtual std::string getDate() = 0;
	virtual std::vector<std::string> getGpuDevices() = 0;
	virtual std::string getKernel() = 0;
	virtual bool getL2NormTestRun() = 0;
	virtual bool getL2NormTestBetweenKernelRun() = 0;
	virtual int getNumberOfTimeSteps() = 0;
	virtual bool getNyTestRun() = 0;
	virtual bool getPhiTestRun() = 0;
	virtual std::vector<double> getResultCheckTime() = 0;
	virtual std::string getSimName() = 0;
	virtual std::vector<int> getSimTime() = 0;

	virtual std::vector<double> getTestTime() = 0;
	virtual std::string getTime() = 0;
	virtual double getViscosity() = 0;
	virtual bool getVTKFileWriting() = 0;

	virtual std::string getFilePath() = 0;
	virtual std::string getSimulationSigniture() = 0;

	virtual std::vector<std::shared_ptr<L2NormLogFileData> > getL2NormLogFileData() = 0;
	virtual std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > getL2NormBetweenKernelsLogFileData() = 0;
	virtual std::vector<std::shared_ptr<NyLogFileData> > getNyLogFileData() = 0;
	virtual std::vector<std::shared_ptr<PhiLogFileData> > getPhiLogFileData() = 0;

	virtual std::shared_ptr<TaylorGreenVortexUxLogFileData> getTaylorGreenVortexUxLogFileData() = 0;
	virtual std::shared_ptr<TaylorGreenVortexUzLogFileData> getTaylorGreenVortexUzLogFileData() = 0;
	virtual std::shared_ptr<ShearWaveLogFileData> getShearWaveLogFileData() = 0;
	
};
#endif