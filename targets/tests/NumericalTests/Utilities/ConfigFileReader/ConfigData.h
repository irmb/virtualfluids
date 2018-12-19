#ifndef CONFIGFILE_DATA_H
#define CONFIGFILE_DATA_H

#include <vector>

struct ConfigDataStruct
{
	std::vector< int> devices;

	std::vector< std::string> kernelsToTest;

	unsigned int numberOfTimeSteps;
	std::vector<double> viscosity;

	std::vector<double> amplitudeTGVux, u0TGVux;
	double l0TGVux;
	std::vector<double> amplitudeTGVuz, v0TGVuz;
	double l0TGVuz;
	std::vector<int> basisTimeStepLengthTGVux, basisTimeStepLengthTGVuz;

	std::vector<double> u0SW, v0SW;
	std::vector<int> basisTimeStepLengthSW;
	double l0SW;

	unsigned int ySliceForCalculation;

	double minOrderOfAccuracy;
	std::vector<std::string> dataToCalcPhiAndNuTest;
	unsigned int startTimeStepCalculationPhiNu;
	unsigned int endTimeStepCalculationPhiNu;
	bool nuAndPhiTest;

	double maxL2NormDiff;
	std::vector<std::string> dataToCalcL2Test;
	unsigned int basicTimeStepL2Norm;
	unsigned int divergentTimeStepL2Norm;
	bool l2NormTest;

	bool l2NormBetweenKernelTest;
	std::string basicKernelL2NormTest;
	std::vector< int> timeStepsL2NormBetweenKernel;
	std::vector< std::string> dataToCalcL2NormBetweenKernel;

	std::vector< bool> tgvUx, tgvUz;
	std::vector< bool> sw;

	unsigned int numberOfGridLevels;
	unsigned int maxLevel;
	std::vector< std::string> grids;

	bool writeFiles;
	std::string filePath;
	unsigned int startStepFileWriter;
	std::string logFilePath;
	
	double rho0;
	std::vector< double> lx;
	std::vector< double> lz;

	bool writeAnalyticalToVTK;

};
#endif 