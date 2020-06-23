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
//! \file CupsAnalyzer.h
//! \ingroup Analyzer
//! \author Stephan Lenz
//=======================================================================================
#ifndef  CupsAnalyzer_H
#define  CupsAnalyzer_H

#include <string>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/Timer/Timer.h"

struct DataBase;

//! \brief Timing analyzer
//! 
//! Measures the update rate of the simulation in Cell Updates Per Second (CUPS) and other timing metrics
class VF_PUBLIC CupsAnalyzer
{
private:
    SPtr<Timer> timer;          //!< Timer for measuring elapse wall clock time
    SPtr<Timer> timerRestart;   //!< Timer for CUPS calculation

    bool outputPerTime;         //!< enables analyzer execution based on elapsed time in seconds

    bool outputPerIter;         //!< enables analyzer execution based on number of iterations

    real outputTime;            //!< controls analyzer execution based on elapsed time   
    uint outputPerTimeCounter;  //!< counts how often the analyzer executed based on the elapsed time

    uint outputIter;            //!< controls analyzer execution based on number of iterations

    unsigned long long numberOfCellUpdatesPerTimeStep;  //!< total number of cell updates that has to be performed for a single time step

    uint counter;               //!< counts the number of executions of this analyzer

public:

    //! constructor
    //! \param dataBase        shared pointer to a \ref DataBase
    //! \param outputPerTime   flag to enable/disable execution based on elapsed time in seconds
    //! \param outputTime      time interval after which the analyzer should be executed
    //! \param outputPerIter   flag to enable/disable execution based on number of iterations
    //! \param outputIter      number of iterations after which the analyzer should be executed
    CupsAnalyzer( SPtr<DataBase> dataBase, 
                  bool outputPerTime = true, real outputTime = 600.0,
                  bool outputPerIter = true, uint outputIter = 10000 );

    //! starts the timing
    void start();

    //! restarts the timing only for CUPS calculation
    void restart();

    //! executes the analyzer
    //! \param iter   current iteration number
    //! \param dt     coarse time step used to compute simulation time
    void run( uint iter, real dt );

private:

    //! \param currentRuntime  amount of elapsed time in seconds
    //! \return  true if the analyzer should be executed
    bool checkOutputPerTime( real currentRuntime );

    //! \param iter current iteration number
    //! \return return true if the analyzer should be executed
    bool checkOutputPerIter( uint iter );

    //! prints the timing information on the logger
    void printCups(uint iter, real simTime, real currentRunTime, real cups);

    //! renders the time into a nice string of format dd:hh:mm:ss.ms
    //! \param time  time in seconds
    //! \return   time as nice string
    std::string getTimeString( real time );
};

#endif
