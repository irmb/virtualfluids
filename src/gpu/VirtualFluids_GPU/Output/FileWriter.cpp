//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ /
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "FileWriter.h"
#include <logger/Logger.h>

#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>

#include <StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"

#include <basics/writer/WbWriterVtkXmlBinary.h>

std::string makePartFileNameEnding(int level, int ID, int part, int timestep)
{
    return "_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(ID) + "_Part_" + StringUtil::toString<int>(part) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";
}

std::string makeCollectionFileNameEnding(int ID, int timestep)
{
    return "_ID_" + StringUtil::toString<int>(ID) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";
}

std::string makePartFileName(const std::string &prefix, int level, int ID, int part, int timestep)
{
    return prefix + "_bin" + makePartFileNameEnding(level, ID, part, timestep);
}

std::string makeMedianPartFileName(const std::string &prefix, int level, int ID, int part, int timestep)
{
    return prefix + "_bin_median" + makePartFileNameEnding(level, ID, part, timestep);
}


std::string makeCollectionFileName(const std::string &prefix, int ID, int timestep)
{
    return prefix + "_bin" + makeCollectionFileNameEnding(ID, timestep);
}

std::string makeMedianCollectionFileName(const std::string &prefix, int ID, int timestep)
{
    return prefix + "_bin_median" + makeCollectionFileNameEnding(ID, timestep);
}

void FileWriter::writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager)
{
    unsigned int timestep = para->getTimestepInit();
    for (int level = para->getCoarse(); level <= para->getFine(); level++) {
        cudaMemoryManager->cudaCopyPrint(level);
        writeTimestep(para, timestep, level);
    }

    this->writeCollectionFile(para, timestep);

    if( para->getCalcMedian() )
        this->writeCollectionFileMedian(para, timestep);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
        writeTimestep(para, timestep, level);

    this->writeCollectionFile(para, timestep);

    if( para->getCalcMedian() )
        this->writeCollectionFileMedian(para, timestep);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level)
{
    const unsigned int numberOfParts = (uint)para->getParH(level)->numberOfNodes / para->getLimitOfNodesForVTK() + 1;
    std::vector<std::string> fnames;
    std::vector<std::string> fnamesMed;

    for (unsigned int i = 1; i <= numberOfParts; i++)
    {
        std::string fname = makePartFileName(para->getFName(), level, para->getMyProcessID(), i, timestep); 
        std::string fnameMed = makeMedianPartFileName(para->getFName(), level, para->getMyProcessID(), i, timestep); 

        fnames.push_back(fname);
        fnamesMed.push_back(fnameMed);
    }

    std::vector<std::string> fnamesLong = writeUnstructuredGridLT(para, level, fnames);
    for(auto fname : fnamesLong)
        this->fileNamesForCollectionFile.push_back(fname.substr( fname.find_last_of('/') + 1 ));


    if (para->getCalcMedian())
    {
        std::vector<std::string> fnamesMedianLong = writeUnstructuredGridMedianLT(para, level, fnamesMed);
        for(auto fname : fnamesMedianLong)
            this->fileNamesForCollectionFileMedian.push_back(fname.substr( fname.find_last_of('/') + 1 ));
    }
}

bool FileWriter::isPeriodicCell(std::shared_ptr<Parameter> para, int level, unsigned int number2, unsigned int number1, unsigned int number3, unsigned int number5)
{
    return (para->getParH(level)->coordinateX[number2] < para->getParH(level)->coordinateX[number1]) ||
           (para->getParH(level)->coordinateY[number3] < para->getParH(level)->coordinateY[number1]) ||
           (para->getParH(level)->coordinateZ[number5] < para->getParH(level)->coordinateZ[number1]);
}

std::vector<std::string> FileWriter::getNodeDataNames(std::shared_ptr<Parameter> para)
{

    std::vector<std::string> nodeDataNames;
    nodeDataNames.push_back("press");
    nodeDataNames.push_back("rho");
    nodeDataNames.push_back("vx1");
    nodeDataNames.push_back("vx2");
    nodeDataNames.push_back("vx3");
    nodeDataNames.push_back("geo");

    if(para->getDiffOn()) 
        nodeDataNames.push_back("conc");

    if(para->getIsBodyForce())
    {
        nodeDataNames.push_back("Fx");
        nodeDataNames.push_back("Fy");
        nodeDataNames.push_back("Fz");
    }

    if(para->getUseTurbulentViscosity())
    {
        nodeDataNames.push_back("nut");
    }

    if (para->getCalcTurbulenceIntensity()) {
        nodeDataNames.push_back("vxx");
        nodeDataNames.push_back("vyy");
        nodeDataNames.push_back("vzz");
        nodeDataNames.push_back("vxy");
        nodeDataNames.push_back("vxz");
        nodeDataNames.push_back("vyz");
    }
    return nodeDataNames;
}

std::vector<std::string> FileWriter::getMedianNodeDataNames(std::shared_ptr<Parameter> para)
{
    std::vector<std::string> nodeDataNames;
    
    if(para->getDiffOn()) 
        nodeDataNames.push_back("conc");
    nodeDataNames.push_back("pressMed");
    nodeDataNames.push_back("rhoMed");
    nodeDataNames.push_back("vx1Med");
    nodeDataNames.push_back("vx2Med");
    nodeDataNames.push_back("vx3Med");
    nodeDataNames.push_back("geo");

    return nodeDataNames;
}

std::string VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFile(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    auto nodeDataNames = this->getNodeDataNames(para);
    std::vector<std::string> cellDataNames;
    std::string pFileName= WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFile, nodeDataNames, cellDataNames);
    this->fileNamesForCollectionFile.clear();
    return pFileName;
}

std::string VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFileMedian(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeMedianCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    std::vector<std::string> nodeDataNames = getMedianNodeDataNames(para);
    std::vector<std::string> cellDataNames;
    std::string pFileName =  WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFileMedian, nodeDataNames, cellDataNames);
    this->fileNamesForCollectionFileMedian.clear();
    return pFileName;
}

std::vector<std::string> FileWriter::writeUnstructuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    std::vector< std::string > nodeDataNames = getNodeDataNames(para);

    std::vector< std::string > outFNames;

    uint dataIndex = 6;

    uint firstConcNode = dataIndex;
    if(para->getDiffOn()) dataIndex++;
    
    uint firstBodyForceNode = dataIndex;
    if(para->getIsBodyForce()) dataIndex+=3;

    uint firstNutNode = dataIndex;
    if(para->getUseTurbulentViscosity()) dataIndex++;

    uint firstTurbulenceNode = dataIndex;
    if (para->getCalcTurbulenceIntensity()) dataIndex += 6;
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    unsigned int startPosition;
    unsigned int endPosition;
    unsigned int sizeOfNodes;
    std::vector< std::vector< double > > nodeData(nodeDataNames.size());

    for (unsigned int part = 0; part < fname.size(); part++)
    {
        if (((part + 1)*para->getLimitOfNodesForVTK()) > (uint)para->getParH(level)->numberOfNodes)
            sizeOfNodes = (uint)para->getParH(level)->numberOfNodes - (part * para->getLimitOfNodesForVTK());
        else
            sizeOfNodes = para->getLimitOfNodesForVTK();

        //////////////////////////////////////////////////////////////////////////
        startPosition = part * para->getLimitOfNodesForVTK();
        endPosition = startPosition + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        for (uint i = 0; i < (uint)nodeDataNames.size(); i++)
            nodeData[i].resize(sizeOfNodes);

        //////////////////////////////////////////////////////////////////////////
        for (unsigned int pos = startPosition; pos < endPosition; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {

                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startPosition;
                neighborsAreFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodeData[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodeData[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodeData[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodeData[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodeData[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodeData[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];

                if(para->getDiffOn())
                    nodeData[firstConcNode][dn1] = (double)para->getParH(level)->concentration[pos];

                if(para->getIsBodyForce())
                {
                    nodeData[firstBodyForceNode    ][dn1] = (double)para->getParH(level)->forceX_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodeData[firstBodyForceNode + 1][dn1] = (double)para->getParH(level)->forceY_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodeData[firstBodyForceNode + 2][dn1] = (double)para->getParH(level)->forceZ_SP[pos] * (double)para->getScaledForceRatio(level);
                }

                if(para->getUseTurbulentViscosity())
                {
                    nodeData[firstNutNode][dn1] = (double)para->getParH(level)->turbViscosity[pos] * (double)para->getScaledViscosityRatio(level);
                }

                if (para->getCalcTurbulenceIntensity()) {
                    nodeData[firstTurbulenceNode    ][dn1] = (double)para->getParH(level)->vxx[pos];
                    nodeData[firstTurbulenceNode + 1][dn1] = (double)para->getParH(level)->vyy[pos];
                    nodeData[firstTurbulenceNode + 2][dn1] = (double)para->getParH(level)->vzz[pos];
                    nodeData[firstTurbulenceNode + 3][dn1] = (double)para->getParH(level)->vxy[pos];
                    nodeData[firstTurbulenceNode + 4][dn1] = (double)para->getParH(level)->vxz[pos];
                    nodeData[firstTurbulenceNode + 5][dn1] = (double)para->getParH(level)->vyz[pos];
                }

                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID)  neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endPosition ||
                    number3 > endPosition ||
                    number4 > endPosition ||
                    number5 > endPosition ||
                    number6 > endPosition ||
                    number7 > endPosition ||
                    number8 > endPosition)  neighborsAreFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startPosition;
                dn3 = number3 - startPosition;
                dn4 = number4 - startPosition;
                dn5 = number5 - startPosition;
                dn6 = number6 - startPosition;
                dn7 = number7 - startPosition;
                dn8 = number8 - startPosition;
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
            }
        }
        outFNames.push_back( WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodeDataNames, nodeData) );
    }
    return outFNames;
}

std::vector<std::string> FileWriter::writeUnstructuredGridMedianLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
    std::vector< std::string > outFNames;

    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    //std::vector< UbTupleUInt8 > cells2;
    std::vector< std::string > nodeDataNames = getMedianNodeDataNames(para);
    int startIndex = para->getDiffOn()? 1 : 0;

    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    unsigned int startPosition;
    unsigned int endPosition;
    unsigned int sizeOfNodes;
    std::vector< std::vector< double > > nodeData(nodeDataNames.size());

    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        if (((part + 1) * para->getLimitOfNodesForVTK()) > (uint)para->getParH(level)->numberOfNodes)
        {
            sizeOfNodes = (uint)para->getParH(level)->numberOfNodes - (part * para->getLimitOfNodesForVTK());
        }
        else
        {
            sizeOfNodes = para->getLimitOfNodesForVTK();
        }
        //////////////////////////////////////////////////////////////////////////
        startPosition = part * para->getLimitOfNodesForVTK();
        endPosition = startPosition + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        for (size_t i = 0; i < nodeDataNames.size(); i++)
            nodeData[i].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos = startPosition; pos < endPosition; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startPosition;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                if(para->getDiffOn())
                    nodeData[0][dn1] = (double)para->getParH(level)->Conc_Med_Out[pos];
                nodeData[startIndex    ][dn1] = para->getParH(level)->press_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodeData[startIndex + 1][dn1] = para->getParH(level)->rho_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodeData[startIndex + 2][dn1] = para->getParH(level)->vx_SP_Med_Out[pos] * para->getVelocityRatio();
                nodeData[startIndex + 3][dn1] = para->getParH(level)->vy_SP_Med_Out[pos] * para->getVelocityRatio();
                nodeData[startIndex + 4][dn1] = para->getParH(level)->vz_SP_Med_Out[pos] * para->getVelocityRatio();
                nodeData[startIndex + 5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
                //////////////////////////////////////////////////////////////////////////
                number2 = para->getParH(level)->neighborX[number1];
                number3 = para->getParH(level)->neighborY[number2];
                number4 = para->getParH(level)->neighborY[number1];
                number5 = para->getParH(level)->neighborZ[number1];
                number6 = para->getParH(level)->neighborZ[number2];
                number7 = para->getParH(level)->neighborZ[number3];
                number8 = para->getParH(level)->neighborZ[number4];
                //////////////////////////////////////////////////////////////////////////
                if (para->getParH(level)->typeOfGridNode[number2] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number3] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number4] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number5] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number6] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number7] != GEO_FLUID ||
                    para->getParH(level)->typeOfGridNode[number8] != GEO_FLUID)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                if (number2 > endPosition ||
                    number3 > endPosition ||
                    number4 > endPosition ||
                    number5 > endPosition ||
                    number6 > endPosition ||
                    number7 > endPosition ||
                    number8 > endPosition)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startPosition;
                dn3 = number3 - startPosition;
                dn4 = number4 - startPosition;
                dn5 = number5 - startPosition;
                dn6 = number6 - startPosition;
                dn7 = number7 - startPosition;
                dn8 = number8 - startPosition;
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid == true) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        outFNames.push_back(WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodeDataNames, nodeData));
        //////////////////////////////////////////////////////////////////////////
    }
    return outFNames;
}