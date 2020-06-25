#ifndef MATHEMATICA_POINT_LIST_IMP_H
#define MATHEMATICA_POINT_LIST_IMP_H

#include "MathematicaPointList.h"

#include <memory>
#include <vector>

class DataPoint;

class MathematicaPointListImp : public MathematicaPointList
{
public:
	static std::shared_ptr<MathematicaPointList> getNewInstance(std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData);

	std::string getListName();

private:
	MathematicaPointListImp();
	MathematicaPointListImp(std::string listName, std::vector<std::shared_ptr<DataPoint> > plotData);

	std::string listName;
};
#endif
