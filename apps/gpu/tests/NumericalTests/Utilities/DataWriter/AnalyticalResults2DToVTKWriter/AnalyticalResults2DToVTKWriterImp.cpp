#include "AnalyticalResults2DToVTKWriterImp.h"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"

#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Utilities/Results/AnalyticalResults/AnalyticalResult.h"
#include <mpi.h>


std::shared_ptr<AnalyticalResults2DToVTKWriterImp> AnalyticalResults2DToVTKWriterImp::getInstance(bool writeAnalyticalResults)
{
	static std::shared_ptr<AnalyticalResults2DToVTKWriterImp> uniqueInstance;
	if (!uniqueInstance)
		uniqueInstance = std::shared_ptr<AnalyticalResults2DToVTKWriterImp>(new AnalyticalResults2DToVTKWriterImp(writeAnalyticalResults));
	return uniqueInstance;
}

AnalyticalResults2DToVTKWriterImp::AnalyticalResults2DToVTKWriterImp(bool writeAnalyticalResults) : writeAnalyticalResults(writeAnalyticalResults)
{

}

void AnalyticalResults2DToVTKWriterImp::writeAnalyticalResult(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult)
{
	if (writeAnalyticalResults) {
		std::cout << "Write Analytical Result To VTK-Files" << std::endl;
		for (int level = para->getCoarse(); level <= para->getFine(); level++) {
#pragma omp parallel for
			for (int timeStep = 0; timeStep < analyticalResult->getNumberOfTimeSteps(); timeStep++) {
				const unsigned int numberOfParts = para->getParH(level)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
				std::vector<std::string> fname;
				unsigned int time = analyticalResult->getTimeSteps().at(timeStep)*analyticalResult->getTimeStepLength();
				for (int j = 1; j <= numberOfParts; j++) {
					std::string filePath = para->getFName();
					filePath.resize(filePath.size() - 5);
					fname.push_back(filePath + "AnalyticalResult/Analytical_cells_bin_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(j) + "_t_" + StringUtil::toString<int>(time) + ".vtk");
				}
				std::cout << "\t Write TimeStep=" << timeStep << " t=" << time << "...";
				writeTimeStep(para, analyticalResult, level, fname, timeStep);
				std::cout << "done." << std::endl;
			}
		}
		std::cout << std::endl;
	}
}


void AnalyticalResults2DToVTKWriterImp::writeTimeStep(std::shared_ptr<Parameter> para, std::shared_ptr<AnalyticalResults> analyticalResult, int level, std::vector<std::string> & fname, int timeStep)
{
	std::vector<UbTupleFloat3 > nodes;
    std::vector<UbTupleUInt8 > cells;
    std::vector<std::string > nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector<std::vector<double > > nodedata(nodedatanames.size());

	maxX = para->getGridX().at(level);
	maxY = para->getGridY().at(level);
	maxZ = para->getGridZ().at(level);

	std::vector<double> press = analyticalResult->getPress()[timeStep];
	std::vector<double> rho = analyticalResult->getRho()[timeStep];
	std::vector<double> vx = analyticalResult->getVx()[timeStep];
	std::vector<double> vy = analyticalResult->getVy()[timeStep];
	std::vector<double> vz = analyticalResult->getVz()[timeStep];

    for (unsigned int part = 0; part < fname.size(); part++){
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

				int numberInResults = CoordResults2DTo1D(x1 - 1.0, x3 - 1.0);
                nodedata[0][dn1] = press[numberInResults];
				nodedata[1][dn1] = rho[numberInResults];
                nodedata[2][dn1] = vx[numberInResults];
                nodedata[3][dn1] = vy[numberInResults];
                nodedata[4][dn1] = vz[numberInResults];
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

int AnalyticalResults2DToVTKWriterImp::CoordResults2DTo1D(int x, int z)
{
	return z * (maxX - 1) + x;
}
