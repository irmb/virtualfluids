#include "FileWriter.h"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include <utilities/StringUtil/StringUtil.h>

#include "Parameter/Parameter.h"

#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <VirtualFluidsBasics/basics/writer/WbWriterVtkXmlBinary.h>

void FileWriter::writeInit(std::shared_ptr<Parameter> para)
{
    unsigned int timestep = para->getTInit();
	for (int level = para->getCoarse(); level <= para->getFine(); level++) {
		para->cudaCopyPrint(level);
		writeTimestep(para, timestep, level);

	}
        
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
        writeTimestep(para, timestep, level);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level)
{
    const unsigned int numberOfParts = para->getParH(level)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
    std::vector<std::string> fname;
    for (unsigned int i = 1; i <= numberOfParts; i++)
        fname.push_back(para->getFName() + "_bin_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(i) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk");
    
    writeUnstrucuredGridLT(para, level, fname);
}

void FileWriter::writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    std::vector< std::string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    double vxmax = 0;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    for (unsigned int part = 0; part < fname.size(); part++)
    {
        if (((part + 1)*para->getlimitOfNodesForVTK()) > para->getParH(level)->size_Mat_SP)
            sizeOfNodes = para->getParH(level)->size_Mat_SP - (part * para->getlimitOfNodesForVTK());
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
            if (para->getParH(level)->geoSP[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordX_SP[pos];
                double x2 = para->getParH(level)->coordY_SP[pos];
                double x3 = para->getParH(level)->coordZ_SP[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsAreFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH(level)->press_SP[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho_SP[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->vx_SP[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->vy_SP[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->vz_SP[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->geoSP[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX_SP[number1];
                number3 = para->getParH(level)->neighborY_SP[number2];
                number4 = para->getParH(level)->neighborY_SP[number1];
                number5 = para->getParH(level)->neighborZ_SP[number1];
                number6 = para->getParH(level)->neighborZ_SP[number2];
                number7 = para->getParH(level)->neighborZ_SP[number3];
                number8 = para->getParH(level)->neighborZ_SP[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->geoSP[number2] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number3] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number4] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number5] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number6] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number7] != GEO_FLUID ||
                    para->getParH(level)->geoSP[number8] != GEO_FLUID)  neighborsAreFluid = false;
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
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
    }
}