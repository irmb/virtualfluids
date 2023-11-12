#ifndef TGV_UX_LOG_FILE_DATA_IMP_H
#define TGV_UX_LOG_FILE_DATA_IMP_H

#include "TaylorGreenVortexUxLogFileData.h"

#include <memory>

class TaylorGreenVortexUxLogFileDataImp : public TaylorGreenVortexUxLogFileData
{
public:
    static std::shared_ptr<TaylorGreenVortexUxLogFileDataImp> getNewInstance();

    std::vector<int> getL0();
    std::vector<double> getUx();
    std::vector<double> getAmplitude();

    void setL0(std::vector<int> l0);
    void setUx(std::vector<double> ux);
    void setAmplitude(std::vector<double> amp);

    ~TaylorGreenVortexUxLogFileDataImp();

private:
    TaylorGreenVortexUxLogFileDataImp();

    std::vector<int> l0;
    std::vector<double> ux;
    std::vector<double> amp;
};
#endif