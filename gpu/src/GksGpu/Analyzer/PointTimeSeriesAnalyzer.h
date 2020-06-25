#ifndef  PointTimeSeriesAnalyzer_H
#define  PointTimeSeriesAnalyzer_H

#include <vector>
#include <string>
#include <memory>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

class  GksMeshAdapter;

namespace GksGpu {

struct DataBase;
struct Parameters;

struct PointTimeSeriesAnalyzerStruct
{
    real* deviceSeries;

    char quantity;

    uint counter;

    uint cellIndex;
};

class VF_PUBLIC PointTimeSeriesAnalyzer
{
public:

    SPtr<DataBase> dataBase;

    uint outputIter;

    real* deviceSeries;

    char quantity;

    uint counter;

    uint cellIndex;

    std::vector<real> hostSeries;

    Vec3 coordinates;

public:

    ~PointTimeSeriesAnalyzer();

    PointTimeSeriesAnalyzer( SPtr<DataBase> dataBase, GksMeshAdapter & adapter, Vec3 coordinates, char quantity, uint outputIter = 10000 );

    void free();

    void allocate();

    void findCellIndex( GksMeshAdapter & adapter, Vec3 coordinates );

    void run( uint iter, Parameters parameters );

    void writeToFile( std::string filename );

    PointTimeSeriesAnalyzerStruct toStruct();

    void download();
};

} // namespace GksGpu

#endif
