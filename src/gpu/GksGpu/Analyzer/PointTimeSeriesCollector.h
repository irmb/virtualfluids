#ifndef  PointTimeSeriesCelloctor_H
#define  PointTimeSeriesCelloctor_H


#include <vector>
#include <string>
#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

class  GksMeshAdapter;

namespace GksGpu {

class  PointTimeSeriesAnalyzer;
struct DataBase;
struct Parameters;

class VIRTUALFLUIDS_GPU_EXPORT PointTimeSeriesCollector
{
public:

    std::vector< SPtr<PointTimeSeriesAnalyzer> > analyzerList;

public:

    ~PointTimeSeriesCollector();

    PointTimeSeriesCollector(  );

    void addAnalyzer( SPtr<DataBase> dataBase, GksMeshAdapter & adapter, Vec3 coordinate, char quantity, uint outputIter = 10000 );

    void run( uint iter, Parameters parameters );

    void writeToFile( std::string filename );
};

} // namespace GksGpu

#endif
