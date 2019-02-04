#ifndef ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_IMP_H
#define ANALYTICAL_RESULTS_2D_TO_VTK_WRITER_IMP_H

#include "AnalyticalResults2DToVTKWriter.h"

#include <vector>


class AnalyticalResults2DToVTKWriterImp : public AnalyticalResults2DToVTKWriter
{
public:
	static std::shared_ptr<AnalyticalResults2DToVTKWriterImp> getInstance(bool writeAnalyticalResults);

	void writeAnalyticalResult(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult);

private:
	AnalyticalResults2DToVTKWriterImp() {};
	AnalyticalResults2DToVTKWriterImp(bool writeAnalyticalResults);

	void writeTimeStep(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult, int level, std::vector<std::string> & fname, int timeStep);
	int CoordResults2DTo1D(int x, int z);

	std::shared_ptr<Parameter> para;
	int maxX, maxY, maxZ;
	bool writeAnalyticalResults;
};
#endif