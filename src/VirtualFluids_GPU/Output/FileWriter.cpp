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
//! \file FileWriter.cpp
//! \ingroup Output
//! \author Martin Schoenherr
//=======================================================================================
#include "FileWriter.h"
#include "GPU/CudaMemoryManager.h"
#include "Parameter/Parameter.h"
#include "VirtualFluidsBasics/basics/writer/WbWriterVtkXmlBinary.h"

void FileWriter::writeInit(SPtr<Parameter> para, SPtr<CudaMemoryManager> cudaManager)
{
    uint timestep = 0;
	uint level = 0;
	cudaManager->cudaCopyDataToHost();
	writeTimestep(para, timestep);
}

void FileWriter::writeTimestep(SPtr<Parameter> para, uint timestep)
{
	const unsigned int numberOfParts = para->getParH()->numberOfNodes / para->getlimitOfNodesForVTK() + 1;
	std::vector<std::string> fname;

	std::ostringstream strTimestep;
	strTimestep << timestep;

	for (unsigned int i = 1; i <= numberOfParts; i++)
	{
		std::ostringstream strPart;
		strPart << i;
		fname.push_back(para->getPathAndFilename() + "_bin_Part_" + strPart.str() + "_t_" + strTimestep.str() + ".vtk");
	}
	writeUnstrucuredGridLT(para, fname);
}

bool FileWriter::isPeriodicCell(SPtr<Parameter> para, uint number2, uint number1, uint number3, uint number5)
{
	return (para->getParH()->coordinateX[number2] < para->getParH()->coordinateX[number1]) ||
		   (para->getParH()->coordinateY[number3] < para->getParH()->coordinateY[number1]) ||
		   (para->getParH()->coordinateZ[number5] < para->getParH()->coordinateZ[number1]);
}

void FileWriter::writeUnstrucuredGridLT(SPtr<Parameter> para, std::vector<std::string >& fname)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleInt8 > cells;
    std::vector< std::string > nodedatanames;
    nodedatanames.push_back("Press");
    nodedatanames.push_back("DRho");
    nodedatanames.push_back("Vx1");
    nodedatanames.push_back("Vx2");
    nodedatanames.push_back("Vx3");
    nodedatanames.push_back("geometry");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    for (unsigned int part = 0; part < fname.size(); part++)
    {
        if (((part + 1)*para->getlimitOfNodesForVTK()) > para->getParH()->numberOfNodes)
            sizeOfNodes = para->getParH()->numberOfNodes - (part * para->getlimitOfNodesForVTK());
        else
            sizeOfNodes = para->getlimitOfNodesForVTK();

        //////////////////////////////////////////////////////////////////////////
        startpos = part * para->getlimitOfNodesForVTK();
        endpos = startpos + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        nodedata[0].resize(sizeOfNodes);
        nodedata[1].resize(sizeOfNodes);
        nodedata[2].resize(sizeOfNodes);
        nodedata[3].resize(sizeOfNodes);
        nodedata[4].resize(sizeOfNodes);
        nodedata[5].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH()->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH()->coordinateX[pos];
                double x2 = para->getParH()->coordinateY[pos];
                double x3 = para->getParH()->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsAreFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH()->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH()->rho[pos] * (double)para->getDensityRatio();
                nodedata[2][dn1] = (double)para->getParH()->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH()->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH()->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH()->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH()->neighborX[number1];
                number3 = para->getParH()->neighborY[number2];
                number4 = para->getParH()->neighborY[number1];
                number5 = para->getParH()->neighborZ[number1];
                number6 = para->getParH()->neighborZ[number2];
                number7 = para->getParH()->neighborZ[number3];
                number8 = para->getParH()->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH()->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH()->typeOfGridNode[number8] != GEO_FLUID)  neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
				//////////////////////////////////////////////////////////////////////////
				if (isPeriodicCell(para, number2, number1, number3, number5))
					continue;
				//////////////////////////////////////////////////////////////////////////
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
    }
}






