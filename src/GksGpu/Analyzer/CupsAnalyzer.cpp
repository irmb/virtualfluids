#include "CupsAnalyzer.h"

#include <cmath>
#include <sstream>
#include <iomanip>

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

        this->printCups( iter, currentRuntime, CUPS );
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

void CupsAnalyzer::printCups(uint iter, real currentRunTime, real cups)
{
    std::stringstream header;
    std::stringstream body;

    header << "| ";
    header << "      Iter" << " | "; 
    header << " runtime/s" << " | "; 
    header << "     MCUPS" << " | ";

    body   << "| ";
    body   << std::setw(10) << std::setprecision(4) << iter  << " | ";
    body   << std::setw(10) << std::setprecision(4) << currentRunTime << " | ";
    body   << std::setw(10) << std::setprecision(4) << cups / 1.0e6 << " | ";

    *logging::out << logging::Logger::INFO_HIGH << "Performance:" << "\n";
    *logging::out << logging::Logger::INFO_HIGH << header.str() << "\n";
    *logging::out << logging::Logger::INFO_HIGH << body.str()   << "\n";
}
