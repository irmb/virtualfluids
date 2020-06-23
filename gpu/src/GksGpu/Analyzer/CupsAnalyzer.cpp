//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __         
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |        
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |        
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |        
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____    
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|   
//      \    \  |    |   ________________________________________________________________    
//       \    \ |    |  |  ______________________________________________________________|   
//        \    \|    |  |  |         __          __     __     __     ______      _______    
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)   
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______    
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  \   
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/   
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can 
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of 
//  the License, or (at your option) any later version.
//  
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT 
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
//  for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file CupsAnalyzer.cpp
//! \ingroup Analyzer
//! \author Stephan Lenz
//=======================================================================================
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
      outputPerTimeCounter(1),
      counter(0)
{
    this->timer        = Timer::makeStart();
    this->timerRestart = Timer::makeStart();

    this->numberOfCellUpdatesPerTimeStep = 0;

    for( uint level = 0; level < dataBase->numberOfLevels; level++ )
    {
        numberOfCellUpdatesPerTimeStep += std::pow( 2, level ) * dataBase->perLevelCount[level].numberOfBulkCells;
    }
}

void CupsAnalyzer::start()
{
    this->counter = 0;
    this->timer->start();
    this->timerRestart->start();
}

void CupsAnalyzer::restart()
{
    this->counter = 0;
    this->timerRestart->start();
}

void CupsAnalyzer::run( uint iter, real dt )
{
    real currentRuntime             = this->timer->getCurrentRuntimeInSeconds();
    real currentRuntimeSinceRestart = this->timerRestart->getCurrentRuntimeInSeconds();

    this->counter++;

    if( checkOutputPerTime(currentRuntime) || checkOutputPerIter(iter) )
    {
        unsigned long long numberOfCellUpdates = this->numberOfCellUpdatesPerTimeStep * (unsigned long long)counter;

        real CUPS = real(numberOfCellUpdates) / currentRuntimeSinceRestart;

        this->printCups( iter, iter * dt, currentRuntime, CUPS );

        this->restart();
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

void CupsAnalyzer::printCups(uint iter, real simTime, real currentRunTime, real cups)
{
    std::stringstream header;
    std::stringstream body;

    header << "| ";
    header << "           Iter" << " | "; 
    header << "      sim. time" << " | "; 
    header << "      wall time" << " | "; 
    header << "          MCUPS" << " | ";

    body   << "| ";
    body   << std::setw(15) << std::setprecision(4) << iter                                        << " | ";
    body   << std::setw(15) << std::setprecision(4) << this->getTimeString(simTime).c_str()        << " | ";
    body   << std::setw(15) << std::setprecision(4) << this->getTimeString(currentRunTime).c_str() << " | ";
    body   << std::setw(15) << std::setprecision(4) << cups / 1.0e6                                << " | ";

    *logging::out << logging::Logger::INFO_HIGH << "Performance:" << "\n";
    *logging::out << logging::Logger::INFO_HIGH << header.str() << "\n";
    *logging::out << logging::Logger::INFO_HIGH << body.str()   << "\n";
}

std::string CupsAnalyzer::getTimeString(real time)
{
    int seconds = int(time);
    int minutes = seconds / 60;
    int hours   = minutes / 60;
    int days    = hours   / 24;

    int milliseconds = int( 1000.0 * ( time - real(seconds)) );

    hours   -=     days * 24;
    minutes -=   ( days * 24 + hours ) * 60;
    seconds -= ( ( days * 24 + hours ) * 60 + minutes ) * 60;

    std::stringstream timeString;
    timeString << std::setw(2) << std::setfill('0') << days    << ":";
    timeString << std::setw(2) << std::setfill('0') << hours   << ":";
    timeString << std::setw(2) << std::setfill('0') << minutes << ":";
    timeString << std::setw(2) << std::setfill('0') << seconds << ".";
    timeString << std::setw(3) << std::setfill('0') << milliseconds;

    return timeString.str();
}


