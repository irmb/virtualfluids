#ifndef ICell_H
#define ICell_H

#include "LBMSystem.h"
#include <vector>

struct ICell3D {
    ICell3D(int size);

    std::vector<LBMReal> TSW;
    std::vector<LBMReal> TNW;
    std::vector<LBMReal> TNE;
    std::vector<LBMReal> TSE;
    std::vector<LBMReal> BSW;
    std::vector<LBMReal> BNW;
    std::vector<LBMReal> BNE;
    std::vector<LBMReal> BSE;
};

inline ICell3D::ICell3D(int size)
{
    TSW.resize(size);
    TNW.resize(size);
    TNE.resize(size);
    TSE.resize(size);
    BSW.resize(size);
    BNW.resize(size);
    BNE.resize(size);
    BSE.resize(size);
}

#endif
