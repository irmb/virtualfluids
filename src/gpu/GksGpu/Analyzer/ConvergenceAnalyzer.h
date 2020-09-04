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
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
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
//! \file ConvergenceAnalyzer.h
//! \ingroup Analyzer
//! \author Stephan Lenz
//=======================================================================================
#ifndef  ConvergenceAnalyzer_H
#define  ConvergenceAnalyzer_H

#include <vector>

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/Timer/Timer.h"

#include "FlowStateData/FlowStateData.cuh"

#include "GksGpu_export.h"

struct DataBase;

//! \brief Analyzes if a simulation is converged and prints the residual change to the logger
class GKSGPU_EXPORT ConvergenceAnalyzer
{
private:

    SPtr<DataBase> dataBase;        //!< shared pointer to the \ref DataBase, which is used to access the simulation data
    
    std::vector<real> dataHostOld;  //!< host copy of the old data
    std::vector<real> dataHostNew;  //!< host copy of the new data

    uint outputIter;                //!< controls analyzer execution based on number of iterations
    
    ConservedVariables convergenceThreshold; //!< if the change is below these thresholds, the run() method returns true

public:

    //! \brief constructor
    //! \param dataBase              shared pointer to a \ref DataBase
    //! \param outputIter            execution interval
    //! \param convergenceThreshold  threshold for considering the simulation as converged
    ConvergenceAnalyzer( SPtr<DataBase> dataBase, uint outputIter = 10000, real convergenceThreshold = 1.0e-6 );

    //! sets the same threshold for all variables
    void setConvergenceThreshold( real convergenceThreshold );

    //! sets different thresholds for different variables
    void setConvergenceThreshold( ConservedVariables convergenceThreshold );

    //! executes the analyzer
    //! \param iter current iteration number, used to decide, whether the analyzer should be executed or not
    //! \return true if simulation is converged, else false
    bool run( uint iter );

private:

    //! writes the formatted residual change to the logger
    //! \param L2Change residual change
    void printL2Change( ConservedVariables L2Change );
};

#endif
