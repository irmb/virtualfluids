#ifndef SHEAR_WAVE_LOG_FILE_DATA_H
#define SHEAR_WAVE_LOG_FILE_DATA_H

#include <vector>

class ShearWaveLogFileData
{
public:
    virtual std::vector<int> getL0() = 0;
    virtual std::vector<double> getUx() = 0;
    virtual std::vector<double> getUz() = 0;
};
#endif