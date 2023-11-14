#ifndef TGV_UZ_LOG_FILE_DATA_IMP_H
#define TGV_UZ_LOG_FILE_DATA_IMP_H

#include "TaylorGreenVortexUzLogFileData.h"

#include <memory>

class TaylorGreenVortexUzLogFileDataImp : public TaylorGreenVortexUzLogFileData
{
public:
    static std::shared_ptr<TaylorGreenVortexUzLogFileDataImp> getNewInstance();

    std::vector<int> getL0();
    std::vector<double> getUz();
    std::vector<double> getAmplitude();

    void setL0(std::vector<int> l0);
    void setUz(std::vector<double> ux);
    void setAmplitude(std::vector<double> amp);

    ~TaylorGreenVortexUzLogFileDataImp();

private:
    TaylorGreenVortexUzLogFileDataImp();

    std::vector<int> l0;
    std::vector<double> ux;
    std::vector<double> amp;
};
#endif