#include "MathematicaPointList.h"

#include "Utilities/DataPoint/DataPoint.h"

std::shared_ptr<MathematicaPointList> MathematicaPointList::getNewInstance(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData)
{
	return std::shared_ptr<MathematicaPointList>(new MathematicaPointList(listName, plotData));
}

std::string MathematicaPointList::getListName()
{
	return listName;
}

MathematicaPointList::MathematicaPointList(std::string listName, std::vector< std::shared_ptr< DataPoint>> plotData) : listName(listName)
{
	mathematicaFunction << listName << "= {";
	for (int i = 0; i < plotData.size(); i++) {
		if (i > 0)
			mathematicaFunction << ", ";
		mathematicaFunction << "{" << plotData.at(i)->getX() << ", " << plotData.at(i)->getY() << "}";
	}

	mathematicaFunction << "}";
}

MathematicaPointList::MathematicaPointList()
{
}


