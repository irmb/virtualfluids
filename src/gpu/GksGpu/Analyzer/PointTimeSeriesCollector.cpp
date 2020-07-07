#include "PointTimeSeriesCollector.h"

#include <iomanip>
#include <fstream>

#include "Core/Logger/Logger.h"

#include "Analyzer/PointTimeSeriesAnalyzer.h"

#include "Parameters/Parameters.h"

namespace GksGpu {

PointTimeSeriesCollector::~PointTimeSeriesCollector()
{
}

PointTimeSeriesCollector::PointTimeSeriesCollector()
{
}

void PointTimeSeriesCollector::addAnalyzer(SPtr<DataBase> dataBase, GksMeshAdapter & adapter, Vec3 coordinate, char quantity, uint outputIter)
{
    auto pointTimeSeriesAnalyzer = std::make_shared<PointTimeSeriesAnalyzer>( dataBase, adapter, coordinate, quantity, outputIter );

    this->analyzerList.push_back( pointTimeSeriesAnalyzer );
}

void PointTimeSeriesCollector::run(uint iter, Parameters parameters)
{
    for( auto analyzer : this->analyzerList )
        analyzer->run(iter, parameters);
}

void PointTimeSeriesCollector::writeToFile(std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "PointTimeSeriesCollector::writeToFile( " << filename << " )" << "\n";

    if( this->analyzerList.size() == 0 )
    {
        *logging::out << logging::Logger::WARNING << "empty!\n";
        return;
    }

    std::ofstream file;

    file.open(filename + ".dat" );

    //////////////////////////////////////////////////////////////////////////

    file << "Number of Points = " << this->analyzerList.size() << "\n";

    for( uint j = 0; j < this->analyzerList.size(); j++ )
    {
        file << "Point " << j << ", ";
        file << "Quantity = "     << this->analyzerList[j]->quantity << ", ";
        file << "Coordinates = ( " << this->analyzerList[j]->coordinates.x << ", "
                                   << this->analyzerList[j]->coordinates.y << ", "
                                   << this->analyzerList[j]->coordinates.z << " )";
        file << "\n";
    }

    //////////////////////////////////////////////////////////////////////////

    uint numberOfTimeSteps = this->analyzerList[0]->hostSeries.size();

    for( uint i = 0; i < numberOfTimeSteps; i++ )
    {
        for( uint j = 0; j < this->analyzerList.size(); j++ )
        {
            file << std::setprecision(15) << this->analyzerList[j]->hostSeries[i] << ", ";
        }

        file << "\n";
    }

    //////////////////////////////////////////////////////////////////////////

    file.close();

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}

} // namespace GksGpu
