#ifndef CONFIGFILE_DATA_H
#define CONFIGFILE_DATA_H

#include "LBM\LB.h"

#include <vector>

struct ConfigDataStruct
{
	std::vector<double> u0SW, v0SW;
	std::vector<double> amplitudeTGV, u0TGV;
	bool nuAndPhiTestTGV, nuAndPhiTestSW;
	bool l2NormTestTGV, l2NormTestSW;
	std::string dataToCalcPhiAndNuTest, dataToCalcL2Test;
	std::vector<double> viscosity;
	real rho0;
	real l0;
	double minOrderOfAccuracy;
	double maxL2NormDiff;
	unsigned int numberOfTimeSteps, basisTimeStepLength;
	unsigned int startStepCalculationPhiNu, endStepCalculationPhiNu;
	unsigned int basicTimeStepL2Norm, divergentTimeStepL2Norm;
	unsigned int startStepFileWriter;
	unsigned int ySliceForCalculation;
	unsigned int maxLevel;
	unsigned int numberOfGridLevels;
	bool writeFiles;
	std::string filePath;
	std::string logFilePath;
	std::vector< std::string> kernelsToTest;
	std::vector< std::string> grids;
	std::vector< real> lx;
	std::vector< real> lz;
	std::vector< int> devices;
	std::vector< bool> tgv;
	std::vector< bool> sw;
};
#endif 