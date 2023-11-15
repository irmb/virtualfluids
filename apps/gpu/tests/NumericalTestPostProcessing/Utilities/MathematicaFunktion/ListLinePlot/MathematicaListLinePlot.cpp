#include "MathematicaListLinePlot.h"

std::shared_ptr<MathematicaListLinePlot> MathematicaListLinePlot::getNewInstance(std::vector< std::shared_ptr< DataPoint>> plotData)
{
    return std::shared_ptr<MathematicaListLinePlot>(new MathematicaListLinePlot(plotData));
}

MathematicaListLinePlot::MathematicaListLinePlot(std::vector< std::shared_ptr< DataPoint>> plotData)
{
    mathematicaFunction << "ListLinePlot[";


    mathematicaFunction << "]";
}

MathematicaListLinePlot::MathematicaListLinePlot()
{
}


