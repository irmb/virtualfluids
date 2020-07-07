#ifndef TGV_UX_LOG_FILE_DATA_H
#define TGV_UX_LOG_FILE_DATA_H

#include <vector>

class TaylorGreenVortexUxLogFileData
{
public:
	virtual std::vector<int> getL0() = 0;
	virtual std::vector<double> getUx() = 0;
	virtual std::vector<double> getAmplitude() = 0;
};
#endif