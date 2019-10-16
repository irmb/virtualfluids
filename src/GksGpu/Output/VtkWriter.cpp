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
//! \file VtkWriter.cpp
//! \ingroup Output
//! \author Stephan Lenz
//=======================================================================================
#include "VtkWriter.h"

#include <vector>
#include <memory>

#include "Core/Logger/Logger.h"

#include "VirtualFluidsBasics/basics/utilities/UbTuple.h"
#include "VirtualFluidsBasics/basics/writer/WbWriterVtkXmlBinary.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

#include "FlowStateData/FlowStateData.cuh"
#include "FlowStateData/FlowStateDataConversion.cuh"
#include "FlowStateData/AccessDeviceData.cuh"

void VtkWriter::write(std::shared_ptr<DataBase> dataBase, Parameters parameters, std::string filename)
{
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << filename << ".vtu" << " ... \n";

    //////////////////////////////////////////////////////////////////////////

    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleInt8  > cells;

    nodes.resize( dataBase->numberOfNodes );
    cells.resize( dataBase->numberOfCells );

    for( uint nodeIdx = 0; nodeIdx < dataBase->numberOfNodes; nodeIdx++ )
    {
        Vec3& node = dataBase->nodeCoordinates[ nodeIdx ];

        nodes[nodeIdx] = makeUbTuple( node.x, node.y, node.z );
    }
    
    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ )
    {
        cells[cellIdx] = makeUbTuple( (int)dataBase->cellToNode[ cellIdx ][ 0 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 1 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 2 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 3 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 4 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 5 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 6 ],
                                      (int)dataBase->cellToNode[ cellIdx ][ 7 ] );
    }

    //////////////////////////////////////////////////////////////////////////

    std::vector< std::string > cellDataNames;
    cellDataNames.push_back("Press");       // 0
    cellDataNames.push_back("Rho");         // 1
    cellDataNames.push_back("Vx");          // 2
    cellDataNames.push_back("Vy");          // 3
    cellDataNames.push_back("Vz");          // 4
    cellDataNames.push_back("Temperature"); // 5
    cellDataNames.push_back("Geometry");    // 6
#ifdef USE_PASSIVE_SCALAR
    cellDataNames.push_back("S_1");         // 7
    cellDataNames.push_back("S_2");         // 8
#endif

    //////////////////////////////////////////////////////////////////////////

    std::vector< std::vector< double > > cellData(cellDataNames.size());

    for( auto& i : cellData ) i.resize( dataBase->numberOfCells );

    for( uint cellIdx = 0; cellIdx < dataBase->numberOfCells; cellIdx++ )
    {
        ConservedVariables cons;

        cons.rho  = dataBase->dataHost[ RHO__(cellIdx, dataBase->numberOfCells) ];
        cons.rhoU = dataBase->dataHost[ RHO_U(cellIdx, dataBase->numberOfCells) ];
        cons.rhoV = dataBase->dataHost[ RHO_V(cellIdx, dataBase->numberOfCells) ];
        cons.rhoW = dataBase->dataHost[ RHO_W(cellIdx, dataBase->numberOfCells) ];
        cons.rhoE = dataBase->dataHost[ RHO_E(cellIdx, dataBase->numberOfCells) ];
#ifdef USE_PASSIVE_SCALAR
        cons.rhoS_1 = dataBase->dataHost[ RHO_S_1(cellIdx, dataBase->numberOfCells) ];
        cons.rhoS_2 = dataBase->dataHost[ RHO_S_2(cellIdx, dataBase->numberOfCells) ];
#endif // USE_PASSIVE_SCALAR

        PrimitiveVariables prim = toPrimitiveVariables(cons, parameters.K);

        real p = 0.5 * prim.rho / prim.lambda;

#ifdef USE_PASSIVE_SCALAR
        real T = getT(prim);
#else // USE_PASSIVE_SCALAR
        real T = 1.0 / prim.lambda;
#endif // USE_PASSIVE_SCALAR

        cellData[0][cellIdx] = p;
        cellData[1][cellIdx] = prim.rho;
        cellData[2][cellIdx] = prim.U;
        cellData[3][cellIdx] = prim.V;
        cellData[4][cellIdx] = prim.W;
        cellData[5][cellIdx] = T;
        cellData[6][cellIdx] = dataBase->isGhostCell(cellIdx);
#ifdef USE_PASSIVE_SCALAR
        cellData[7][cellIdx] = prim.S_1;
        cellData[8][cellIdx] = prim.S_2;
#endif
    }

    //////////////////////////////////////////////////////////////////////////

    WbWriterVtkXmlBinary::getInstance()->writeOctsWithCellData(filename, nodes, cells, cellDataNames, cellData);

    //////////////////////////////////////////////////////////////////////////

    *logging::out << logging::Logger::INFO_INTERMEDIATE << "done!\n";
}
