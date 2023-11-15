#ifndef ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_H
#define ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_H

#include <memory>
#include <string>

class Parameter;
class AnalyticalResults;

class AnalyticalResults2DToVTKWriter
{
public:
    virtual ~AnalyticalResults2DToVTKWriter() = default;
    virtual void writeAnalyticalResult(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult) = 0;

};
#endif