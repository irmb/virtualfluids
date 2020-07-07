#ifndef SHEAR_WAVE_LOG_FILE_DATA_IMP_H
#define SHEAR_WAVE_LOG_FILE_DATA_IMP_H

#include "ShearWaveLogFileData.h"

#include <memory>

class ShearWaveLogFileDataImp : public ShearWaveLogFileData
{
public:
	static std::shared_ptr<ShearWaveLogFileDataImp> getNewInstance();

	std::vector<int> getL0();
	std::vector<double> getUx();
	std::vector<double> getUz();

	void setL0(std::vector<int> l0);
	void setUx(std::vector<double> ux);
	void setUz(std::vector<double> amp);

	~ShearWaveLogFileDataImp();

private:
	ShearWaveLogFileDataImp();

	std::vector<int> l0;
	std::vector<double> ux;
	std::vector<double> uz;
};
#endif