#ifndef TGV_UZ_LOG_FILE_DATA_H
#define TGV_UZ_LOG_FILE_DATA_H

#include <vector>

class TaylorGreenVortexUzLogFileData
{
public:
	virtual std::vector<int> getL0() = 0;
	virtual std::vector<double> getUz() = 0;
	virtual std::vector<double> getAmplitude() = 0;
};
#endif