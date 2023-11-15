#include "MathematicaListPlotImp.h"

#include "Utilities/MathematicaFunction/PointList/MathematicaPointList.h"

std::shared_ptr<MathematicaListPlot> MathematicaListPlotImp::getNewInstance(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
    return std::shared_ptr<MathematicaListPlotImp>(new MathematicaListPlotImp(pointList, plotType, xAxes, yAxes));
}

MathematicaListPlotImp::MathematicaListPlotImp(std::vector<std::shared_ptr<MathematicaPointList>> pointList, std::string plotType, std::string xAxes, std::string yAxes)
{
    mathematicaFunction << plotType << "[{";
    for (int i = 0; i < pointList.size(); i++) {
        if(i < pointList.size() - 1)
            mathematicaFunction << pointList.at(i)->getListName() << ", ";
        else
            mathematicaFunction << pointList.at(i)->getListName() << "}";
    }
    mathematicaFunction << ", PlotLegends -> {\"";
    for (int i = 0; i < pointList.size(); i++) {
        if (i < pointList.size() - 1)
            mathematicaFunction << pointList.at(i)->getListName() << "\", \"";
        else
            mathematicaFunction << pointList.at(i)->getListName() << "\"}";
    }
    mathematicaFunction << ", AxesLabel -> {\"" << xAxes << "\", \"" << yAxes << "\"}, Joined -> True, PlotMarkers->Automatic, PlotStyle -> Dashed]";
}

MathematicaListPlotImp::MathematicaListPlotImp()
{
}