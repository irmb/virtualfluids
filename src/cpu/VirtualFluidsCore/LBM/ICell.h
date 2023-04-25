#ifndef ICell_H
#define ICell_H

#include "LBMSystem.h"
#include <vector>

struct ICell3D {
    ICell3D(int size);

    std::vector<real> TSW;
    std::vector<real> TNW;
    std::vector<real> TNE;
    std::vector<real> TSE;
    std::vector<real> BSW;
    std::vector<real> BNW;
    std::vector<real> BNE;
    std::vector<real> BSE;
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
