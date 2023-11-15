#include "MathematicaPointListImp.h"

#include "Utilities/DataPoint/DataPoint.h"
#include "MathematicaPointList.h"

#include <iomanip>
#include <limits>

std::shared_ptr<MathematicaPointList> MathematicaPointListImp::getNewInstance(std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData)
{
    return std::shared_ptr<MathematicaPointList>(new MathematicaPointListImp(listName, plotData));
}

std::string MathematicaPointListImp::getListName()
{
    return listName;
}

MathematicaPointListImp::MathematicaPointListImp(std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData) : listName(listName)
{
    mathematicaFunction << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    mathematicaFunction <<listName <<"= {";
    for (int i = 0; i < plotData.size(); i++) {
        if (i > 0)
            mathematicaFunction <<", ";
        mathematicaFunction <<"{" << plotData.at(i)->getX() <<", " << plotData.at(i)->getY() <<"}";
    }

    mathematicaFunction << "};";
}

MathematicaPointListImp::MathematicaPointListImp()
{
}