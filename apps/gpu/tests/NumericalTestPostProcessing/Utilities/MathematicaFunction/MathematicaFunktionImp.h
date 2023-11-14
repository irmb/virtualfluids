#ifndef MATHEMATICA_FUNCTION_IMP_H
#define MATHEMATICA_FUNCTION_IMP_H

#include "MathematicaFunktion.h"
#include <sstream>

class MathematicaFunctionImp : public MathematicaFunction
{
public:
    std::string getFunction();

protected:
    MathematicaFunctionImp();
    std::ostringstream mathematicaFunction;

};
#endif