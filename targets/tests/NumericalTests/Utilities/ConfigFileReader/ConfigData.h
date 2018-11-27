#ifndef CONFIGFILE_DATA_H
#define CONFIGFILE_DATA_H

#include <vector>

struct ConfigDataStruct
{
	std::vector< int> devices;

	std::vector< std::string> kernelsToTest;

	unsigned int numberOfTimeSteps;
	unsigned int basisTimeStepLength;
	std::vector<double> viscosity;

	std::vector<double> amplitudeTGV, u0TGV;
	std::vector<double> u0SW, v0SW;
	
	unsigned int ySliceForCalculation;

	double minOrderOfAccuracy;
	std::vector<std::string> dataToCalcPhiAndNuTest;
	unsigned int startTimeStepCalculationPhiNu;
	unsigned int endTimeStepCalculationPhiNu;
	bool nuAndPhiTestTGV;
	bool nuAndPhiTestSW;

	double maxL2NormDiff;
	std::vector<std::string> dataToCalcL2Test;
	unsigned int basicTimeStepL2Norm;
	unsigned int divergentTimeStepL2Norm;
	bool l2NormTestTGV;
	bool l2NormTestSW;

	std::vector< bool> tgv;
	std::vector< bool> sw;

	unsigned int numberOfGridLevels;
	unsigned int maxLevel;
	std::vector< std::string> grids;

	bool writeFiles;
	std::string filePath;
	unsigned int startStepFileWriter;
	std::string logFilePath;
	
	double rho0;
	double l0;
	std::vector< double> lx;
	std::vector< double> lz;

};
#endif 