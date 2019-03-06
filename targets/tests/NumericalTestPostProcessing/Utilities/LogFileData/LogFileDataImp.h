#ifndef LOGFILE_DATA_IMP_H
#define LOGFILE_DATA_IMP_H

#include "LogFileData.h"

#include <memory>

class LogFileDataImp : public LogFileData
{
public:
	static std::shared_ptr<LogFileDataImp> getNewInstance();

	std::vector<int> getAnalyticalVTKWritingTime();
	std::vector<double> getBasicGridLengths();
	int getBasisTimeStepLength();
	std::string getDate();
	std::vector<std::string> getGpuDevices();
	std::string getKernel();
	bool getL2NormTestRun();
	bool getL2NormTestBetweenKernelRun();
	int getNumberOfTimeSteps();
	bool getNyTestRun();
	bool getPhiTestRun();
	std::vector<double> getResultCheckTime();
	std::string getSimName();
	std::vector<int> getSimTime();
	std::vector<double> getTestTime();
	std::string getTime();
	double getViscosity();
	bool getVTKFileWriting();
	std::string getFilePath();
	std::string getSimulationSigniture();

	std::vector<std::shared_ptr<L2NormLogFileData> > getL2NormLogFileData();
	std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > getL2NormBetweenKernelsLogFileData();
	std::vector<std::shared_ptr<PhiLogFileData> > getPhiLogFileData();
	std::vector<std::shared_ptr<NyLogFileData> > getNyLogFileData();
	
	std::shared_ptr<TaylorGreenVortexUxLogFileData> getTaylorGreenVortexUxLogFileData();
	std::shared_ptr<TaylorGreenVortexUzLogFileData> getTaylorGreenVortexUzLogFileData();
	std::shared_ptr<ShearWaveLogFileData> getShearWaveLogFileData();

	void setAnalyticalVTKWritingTime(std::vector<int> analyticalVTKWritingTime);
	void setBasicGridLengths(std::vector<double> basicGridLenghts);
	void setBasisTimeStepLength(int basisTimeStepLength);
	void setDate(std::string date);
	void setGpuDevices(std::vector<std::string> gpuDevices);
	void setKernel(std::string kernelName);
	void setL2NormTestRun(bool l2NormTestRun);
	void setL2NormTestBetweenKernelRun(bool l2NormTestBetweenKernelRun);
	void setNumberOfTimeSteps(int numberOfTimeSteps);
	void setPhiTestRun(bool phiTestRun);
	void setNyTestRun(bool nyTestRun);
	void setResultCheckTime(std::vector<double> resultsCheckTime);
	void setSimName(std::string simName);
	void setSimTime(std::vector<int> simTime);
	void setTestTime(std::vector<double> testTime);
	void setTime(std::string time);
	void setViscosity(double viscosity);
	void setVTKFileWriting(bool vtkFileWriting);
	void setFilePath(std::string filePath);
	void setSimulationSigniture(std::string simulationSigniture);

	void setL2NormLogFileData(std::vector<std::shared_ptr<L2NormLogFileData> > data);
	void setL2NormBetweenKernelsLogFileData(std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > data);
	void setNyLogFileData(std::vector<std::shared_ptr<NyLogFileData> > data);
	void setPhiLogFileData(std::vector<std::shared_ptr<PhiLogFileData> > data);
	
	void setTaylorGreenVortexUxLogFileData(std::shared_ptr<TaylorGreenVortexUxLogFileData> data);
	void setTaylorGreenVortexUzLogFileData(std::shared_ptr<TaylorGreenVortexUzLogFileData> data);
	void setShearWaveLogFileData(std::shared_ptr<ShearWaveLogFileData> data);

	~LogFileDataImp();

private:
	LogFileDataImp();

	std::vector<int> analyticalVTKWritingTime;
	std::vector<double> basicGridLenghts;
	int basisTimeStepLength;
	std::string date;
	std::vector<std::string> gpuDevices;
	std::string kernelName;
	bool l2NormTestRun;
	bool l2NormTestBetweenKernelRun;
	int numberOfTimeSteps;
	bool phiTestRun;
	bool nyTestRun;
	std::vector<double> resultsCheckTime;
	std::string simName;
	std::vector<int> simTime;
	std::vector<double> testTime;
	std::string time;
	double viscosity;
	bool vtkFileWriting;
	std::string filePath;
	std::string simulationSigniture;

	std::vector<std::shared_ptr<L2NormBetweenKernelsLogFileData> > l2NormBetweenKernelsDataLogFileData;
	std::vector<std::shared_ptr<PhiLogFileData> > phiLogFileData;
	std::vector<std::shared_ptr<NyLogFileData> > nyLogFileData;
	std::vector<std::shared_ptr<L2NormLogFileData> > l2NormLogFileData;

	std::shared_ptr<ShearWaveLogFileData> shearWaveLogFileData;
	std::shared_ptr<TaylorGreenVortexUxLogFileData> tgvUxLogFileData;
	std::shared_ptr<TaylorGreenVortexUzLogFileData> tgvUzLogFileData;
};
#endif