#ifndef MATHEMATICA_LIST_OF_LISTS_H
#define MATHEMATICA_LIST_OF_LISTS_H

#include "../MathematicaFunktionImp.h"

class MathematicaListOfLists : public MathematicaFunctionImp
{
public:
	virtual std::string getListName() = 0;
};
#endif