#ifndef MATHEMATICA_PLOT_H
#define MATHEMATICA_PLOT_H

#include "Utilities/MathematicaFunction/MathematicaFunktionImp.h"

#include <memory>
#include <vector>

class DataPoint;

class MathematicaListLinePlot : public MathematicaFunctionImp
{
public:
    static std::shared_ptr< MathematicaListLinePlot> getNewInstance(std::vector< std::shared_ptr< DataPoint>> plotData);
    
private:
    MathematicaListLinePlot();
    MathematicaListLinePlot(std::vector< std::shared_ptr< DataPoint>> plotData);

};
#endif
