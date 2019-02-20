#ifndef MATHEMATICA_PLOT_IMP_H
#define MATHEMATICA_PLOT_IMP_H

#include "MathematicaListPlot.h"

#include <memory>
#include <vector>

class MathematicaPointList;

class MathematicaListPlotImp : public MathematicaListPlot
{
public:
	static std::shared_ptr<MathematicaListPlot> getNewInstance(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes);

private:
	MathematicaListPlotImp();
	MathematicaListPlotImp(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes);
};
#endif
