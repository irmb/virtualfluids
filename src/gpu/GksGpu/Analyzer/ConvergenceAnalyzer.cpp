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
//! \file ConvergenceAnalyzer.cpp
//! \ingroup Analyzer
//! \author Stephan Lenz
//=======================================================================================
#include "ConvergenceAnalyzer.h"

#include <cmath>
#include <sstream>
#include <iomanip>

#include "Core/Logger/Logger.h"

#include "DataBase/DataBase.h"

#include "FlowStateData/AccessDeviceData.cuh"

ConvergenceAnalyzer::ConvergenceAnalyzer(SPtr<DataBase> dataBase, uint outputIter, real convergenceThreshold)
{
    this->dataBase = dataBase;

    this->outputIter = outputIter;

    this->setConvergenceThreshold( convergenceThreshold );

    this->dataHostOld = dataBase->dataHost;
    this->dataHostNew = dataBase->dataHost;
}

void ConvergenceAnalyzer::setConvergenceThreshold(real convergenceThreshold)
{
    this->convergenceThreshold.rho  = convergenceThreshold;
    this->convergenceThreshold.rhoU = convergenceThreshold;
    this->convergenceThreshold.rhoV = convergenceThreshold;
    this->convergenceThreshold.rhoW = convergenceThreshold;
    this->convergenceThreshold.rhoE = convergenceThreshold;
#ifdef USE_PASSIVE_SCALAR
    this->convergenceThreshold.rhoS_1 = convergenceThreshold;
    this->convergenceThreshold.rhoS_2 = convergenceThreshold;
#endif //USE_PASSIVE_SCALAR
}

void ConvergenceAnalyzer::setConvergenceThreshold(ConservedVariables convergenceThreshold)
{
    this->convergenceThreshold.rho  = convergenceThreshold.rho ;
    this->convergenceThreshold.rhoU = convergenceThreshold.rhoU;
    this->convergenceThreshold.rhoV = convergenceThreshold.rhoV;
    this->convergenceThreshold.rhoW = convergenceThreshold.rhoW;
    this->convergenceThreshold.rhoE = convergenceThreshold.rhoE;
#ifdef USE_PASSIVE_SCALAR
    this->convergenceThreshold.rhoS_1 = convergenceThreshold.rhoS_1;
    this->convergenceThreshold.rhoS_2 = convergenceThreshold.rhoS_2;
#endif //USE_PASSIVE_SCALAR
}

bool ConvergenceAnalyzer::run(uint iter)
{
    if( iter % outputIter != 0 ) return false;

    this->dataBase->copyDataDeviceToHost( this->dataHostNew.data() );

    ConservedVariables changeSquareSum, consSquareSum;

    for( uint cellIdx = 0; cellIdx < this->dataBase->numberOfCells; cellIdx++  ){

        ConservedVariables change, cons;

        cons.rho  = this->dataHostNew[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = this->dataHostNew[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = this->dataHostNew[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = this->dataHostNew[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = this->dataHostNew[ RHO_E(cellIdx, dataBase->numberOfCells) ];

        change.rho  = cons.rho  - this->dataHostOld[ RHO__(cellIdx, dataBase->numberOfCells) ];
        change.rhoU = cons.rhoU - this->dataHostOld[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        change.rhoV = cons.rhoV - this->dataHostOld[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        change.rhoW = cons.rhoW - this->dataHostOld[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        change.rhoE = cons.rhoE - this->dataHostOld[ RHO_E(cellIdx, dataBase->numberOfCells) ];
    
        changeSquareSum.rho  += change.rho  * change.rho ;
        changeSquareSum.rhoU += change.rhoU * change.rhoU;
        changeSquareSum.rhoV += change.rhoV * change.rhoV;
        changeSquareSum.rhoW += change.rhoW * change.rhoW;
        changeSquareSum.rhoE += change.rhoE * change.rhoE;
    
        consSquareSum.rho  += cons.rho  * cons.rho ;
        consSquareSum.rhoU += cons.rhoU * cons.rhoU;
        consSquareSum.rhoV += cons.rhoV * cons.rhoV;
        consSquareSum.rhoW += cons.rhoW * cons.rhoW;
        consSquareSum.rhoE += cons.rhoE * cons.rhoE;
    }

    ConservedVariables L2Change;

    L2Change.rho  = std::sqrt( changeSquareSum.rho  / consSquareSum.rho  );
    L2Change.rhoU = std::sqrt( changeSquareSum.rhoU / consSquareSum.rhoU );
    L2Change.rhoV = std::sqrt( changeSquareSum.rhoV / consSquareSum.rhoV );
    L2Change.rhoW = std::sqrt( changeSquareSum.rhoW / consSquareSum.rhoW );
    L2Change.rhoE = std::sqrt( changeSquareSum.rhoE / consSquareSum.rhoE );

    this->dataHostOld = this->dataHostNew;

    this->printL2Change( L2Change );

    if( L2Change.rho  < this->convergenceThreshold.rho  &&
        L2Change.rhoU < this->convergenceThreshold.rhoU &&
        L2Change.rhoV < this->convergenceThreshold.rhoV &&
        L2Change.rhoW < this->convergenceThreshold.rhoW &&
        L2Change.rhoE < this->convergenceThreshold.rhoE )
    {
        return true;
    }

    return false;
}

void ConvergenceAnalyzer::printL2Change(ConservedVariables L2Change)
{
    std::stringstream header;
    std::stringstream body;

    header << "| ";
    header << "       rho" << " | "; 
    header << "      rhoU" << " | "; 
    header << "      rhoV" << " | "; 
    header << "      rhoW" << " | "; 
    header << "      rhoE" << " | ";

    body   << "| ";
    body   << std::setw(10) << std::setprecision(4) << L2Change.rho  << " | ";
    body   << std::setw(10) << std::setprecision(4) << L2Change.rhoU << " | ";
    body   << std::setw(10) << std::setprecision(4) << L2Change.rhoV << " | ";
    body   << std::setw(10) << std::setprecision(4) << L2Change.rhoW << " | ";
    body   << std::setw(10) << std::setprecision(4) << L2Change.rhoE << " | ";

    *logging::out << logging::Logger::INFO_HIGH << "Residual L2-Change:" << "\n";
    *logging::out << logging::Logger::INFO_HIGH << header.str() << "\n";
    *logging::out << logging::Logger::INFO_HIGH << body.str()   << "\n";
}
