#ifndef L2_NORM_LOG_FILE_DATA_IMP_H
#define L2_NORM_LOG_FILE_DATA_IMP_H

#include "L2NormLogFileData.h"

#include <memory>

class L2NormLogFileDataImp : public L2NormLogFileData
{
public:
	static std::shared_ptr<L2NormLogFileDataImp> getNewInstance();

	std::vector<double> getBasicGridLengths();
	std::string getDataToCalc();
	std::string getNormalizeData();
	int getBasicTimeStep();
	int getDivergentTimeStep();
	std::vector<double> getL2NormForBasicTimeStep();
	std::vector<double> getL2NormForDivergentTimeStep();
	std::vector<double> getL2NormDiff();

	void setBasicGridLengths(std::vector<double> basicGridLengths);
	void setDataToCalc(std::string dataToCalc);
	void setNormalizeData(std::string normalizeData);
	void setBasicTimeStep(int basicTimeStep);
	void setDivergentTimeStep(int divergentTimeStep);
	void setL2NormForBasicTimeStep(std::vector<double> l2NormForBasicTimeStep);
	void setL2NormForDivergentTimeStep(std::vector<double> l2NormForDivergentTimeStep);
	void setL2NormDiff(std::vector<double> l2NormDiff);

	~L2NormLogFileDataImp();

private:
	L2NormLogFileDataImp();

	std::vector<double> basicGridLengths;
	std::string dataToCalc;
	std::string normalizeData;
	int basicTimeStep;
	int divergentTimeStep;
	std::vector<double> l2NormForBasicTimeStep;
	std::vector<double> l2NormForDivergentTimeStep;
	std::vector<double> l2NormDiff;
};
#endif