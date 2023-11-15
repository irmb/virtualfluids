#ifndef MATHEMATICA_POINT_LIST_H
#define MATHEMATICA_POINT_LIST_H

#include "../MathematicaFunktionImp.h"

class MathematicaPointList : public MathematicaFunctionImp
{
public:
    virtual std::string getListName() = 0;
};
#endif
