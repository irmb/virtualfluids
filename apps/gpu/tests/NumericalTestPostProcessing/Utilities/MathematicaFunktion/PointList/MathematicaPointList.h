#ifndef MATHEMATICA_POINT_LIST_H
#define MATHEMATICA_POINT_LIST_H

#include "Utilities/MathematicaFunction/MathematicaFunktionImp.h"

#include <memory>
#include <vector>

class DataPoint;

class MathematicaPointList : public MathematicaFunctionImp
{
public:
    static std::shared_ptr< MathematicaPointList> getNewInstance(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData);

    std::string getListName();

private:
    MathematicaPointList();
    MathematicaPointList(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData);


    std::string listName;
};
#endif
