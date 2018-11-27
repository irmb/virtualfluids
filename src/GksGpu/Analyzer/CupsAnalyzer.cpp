#include "CupsAnalyzer.h"

#include <cmath>

#include "Core/Logger/Logger.h"

#include "VirtualFluidsDefinitions.h"

#include "DataBase/DataBase.h"

CupsAnalyzer::CupsAnalyzer(SPtr<DataBase> dataBase, 
                           bool outputPerTime, real outputTime, 
                           bool outputPerIter, uint outputIter)
    : outputPerTime(outputPerTime),
      outputTime(outputTime),
      outputPerIter(outputPerIter),
      outputIter(outputIter),
      outputPerTimeCounter(1)
{
    this->numberOfCellUpdatesPerTimeStep = 0;

    for( uint level = 0; level < dataBase->numberOfLevels; level++ )
    {
        numberOfCellUpdatesPerTimeStep += std::pow( 2, level ) * dataBase->perLevelCount[level].numberOfBulkCells;
    }
}

void CupsAnalyzer::start()
{
    this->timer.start();
}

void CupsAnalyzer::run( uint iter )
{
    real currentRuntime = this->timer.getCurrentRuntimeInSeconds();

    if( checkOutputPerTime(currentRuntime) || checkOutputPerIter(iter) )
    {
        unsigned long long numberOfCellUpdates = this->numberOfCellUpdatesPerTimeStep * (unsigned long long)iter;

        real CUPS = real(numberOfCellUpdates) / currentRuntime;

        *logging::out << logging::Logger::INFO_HIGH << "Iteration:   " << iter << "\n";
        *logging::out << logging::Logger::INFO_HIGH << "Run time:    " << currentRuntime << " s\n";
        *logging::out << logging::Logger::INFO_HIGH << "Update rate: " << CUPS / 1.0e6 << " MCUPS\n";
    }

    if( checkOutputPerTime(currentRuntime) )
    {
        outputPerTimeCounter++;
    }
}

bool CupsAnalyzer::checkOutputPerTime(real currentRuntime)
{
    return outputPerTime && ( ( currentRuntime - outputPerTimeCounter * outputTime ) > 0 );
}

bool CupsAnalyzer::checkOutputPerIter(uint iter)
{
    return outputPerIter && (iter % outputIter == 0);
}
