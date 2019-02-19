#ifndef L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_IMP_H
#define L2NORM_BETWEEN_KERNELS_LOG_FILE_DATA_IMP_H

#include "L2NormBetweenKernelsLogFileData.h"

#include <memory>

class L2NormBetweenKernelsLogFileDataImp : public L2NormBetweenKernelsLogFileData
{
public:
	static std::shared_ptr<L2NormBetweenKernelsLogFileDataImp> getNewInstance();

	std::vector<double> getBasicGridLengths();
	std::string getBasicKernel();
	std::string getDataToCalculate();
	int getTimeStep();
	std::vector<double> getL2NormForBasicKernel();
	std::vector<double> getL2NormForDivergentKernel();
	std::vector<double> getL2NormBetweenKernels();
	std::string getNormalizeData();

	void setBasicGridLengths(std::vector<double> basicGridLengths);
	void setBasicKernel(std::string basicKernel);
	void setDataToCalculate(std::string dataToCalc);
	void setTimeStep(int timeStep);
	void setL2NormForBasicKernel(std::vector<double> l2Norm);
	void setL2NormForDivergentKernel(std::vector<double> l2Norm);
	void setL2NormBetweenKernels(std::vector<double> l2Norm);
	void setNormalizeData(std::string normalizeData);

	~L2NormBetweenKernelsLogFileDataImp();

private:
	L2NormBetweenKernelsLogFileDataImp();

	std::vector<double> basicGridLengths;
	std::string basicKernel;
	std::string dataToCalc;
	int timeStep;
	std::vector<double> l2NormForBasicKernel;
	std::vector<double> l2NormForDivergentKernel;
	std::vector<double> l2NormBetweenKernels;
	std::string normalizeData;
};
#endif