//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ /
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "FileWriter.h"
#include <logger/Logger.h>

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include <Core/StringUtilities/StringUtil.h>

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

std::string makePartFileName(std::string prefix, int level, int ID, int part, int timestep)
{
    return prefix + "_bin" + makePartFileNameEnding(level, ID, part, timestep);
}

std::string makeMedianPartFileName(std::string prefix, int level, int ID, int part, int timestep)
{
    return prefix + "_bin_median" + makePartFileNameEnding(level, ID, part, timestep);
}


std::string makeCollectionFileName(std::string prefix, int ID, int timestep)
{
    return prefix + "_bin" + makeCollectionFileNameEnding(ID, timestep);
}

std::string makeMedianCollectionFileName(std::string prefix, int ID, int timestep)
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
    const unsigned int numberOfParts = (uint)para->getParH(level)->numberOfNodes / para->getlimitOfNodesForVTK() + 1;
    std::vector<std::string> fnames;
    std::vector<std::string> fnamesMed;

    for (unsigned int i = 1; i <= numberOfParts; i++)
    {
        std::string fname = makePartFileName(para->getFName(), level, para->getMyProcessID(), i, timestep); 
        std::string fnameMed = makeMedianPartFileName(para->getFName(), level, para->getMyProcessID(), i, timestep); 

        fnames.push_back(fname);
        fnamesMed.push_back(fnameMed);
    }

    std::vector<std::string> fnamesLong = writeUnstrucuredGridLT(para, level, fnames);
    for(auto fname : fnamesLong)
        this->fileNamesForCollectionFile.push_back(fname.substr( fname.find_last_of('/') + 1 ));


    if (para->getCalcMedian())
    {
        std::vector<std::string> fnamesMedianLong = writeUnstrucuredGridMedianLT(para, level, fnamesMed);
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

    std::vector<std::string> nodedatanames;
    nodedatanames.push_back("press");
    nodedatanames.push_back("rho");
    nodedatanames.push_back("vx1");
    nodedatanames.push_back("vx2");
    nodedatanames.push_back("vx3");
    nodedatanames.push_back("geo");

    if(para->getDiffOn()) 
        nodedatanames.push_back("conc");

    if(para->getIsBodyForce())
    {
        nodedatanames.push_back("Fx");
        nodedatanames.push_back("Fy");
        nodedatanames.push_back("Fz");
    }

    if(para->getUseTurbulentViscosity())
    {
        nodedatanames.push_back("nut");
    }

    if (para->getCalcTurbulenceIntensity()) {
        nodedatanames.push_back("vxx");
        nodedatanames.push_back("vyy");
        nodedatanames.push_back("vzz");
        nodedatanames.push_back("vxy");
        nodedatanames.push_back("vxz");
        nodedatanames.push_back("vyz");
    }
    return nodedatanames;
}

std::vector<std::string> FileWriter::getMedianNodeDataNames(std::shared_ptr<Parameter> para)
{
    std::vector<std::string> nodedatanames;
    
    if(para->getDiffOn()) 
        nodedatanames.push_back("conc");
    nodedatanames.push_back("pressMed");
    nodedatanames.push_back("rhoMed");
    nodedatanames.push_back("vx1Med");
    nodedatanames.push_back("vx2Med");
    nodedatanames.push_back("vx3Med");
    nodedatanames.push_back("geo");

    return nodedatanames;
}

std::string VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFile(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    auto nodedatanames = this->getNodeDataNames(para);
    std::vector<std::string> celldatanames;
    std::string pFileName= WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFile, nodedatanames, celldatanames);
    this->fileNamesForCollectionFile.clear();
    return pFileName;
}

std::string VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFileMedian(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeMedianCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    std::vector<std::string> nodedatanames = getMedianNodeDataNames(para);
    std::vector<std::string> celldatanames;
    std::string pFileName =  WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFileMedian, nodedatanames, celldatanames);
    this->fileNamesForCollectionFileMedian.clear();
    return pFileName;
}

std::vector<std::string> FileWriter::writeUnstrucuredGridLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    std::vector< std::string > nodedatanames = getNodeDataNames(para);

    std::vector< std::string > outFNames;

    uint dataIndex = 6;

    uint firstConcNode = dataIndex;
    if(para->getDiffOn()) dataIndex++;
    
    uint firstBodyForceNode = dataIndex;
    if(para->getIsBodyForce()) dataIndex+=3;

    uint firstNutNode = dataIndex;
    if(para->getUseTurbulentViscosity()) dataIndex++;

    uint firstTurbNode = dataIndex;
    if (para->getCalcTurbulenceIntensity()) dataIndex += 6;
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    for (unsigned int part = 0; part < fname.size(); part++)
    {
        if (((part + 1)*para->getlimitOfNodesForVTK()) > (uint)para->getParH(level)->numberOfNodes)
            sizeOfNodes = (uint)para->getParH(level)->numberOfNodes - (part * para->getlimitOfNodesForVTK());
        else
            sizeOfNodes = para->getlimitOfNodesForVTK();

        //////////////////////////////////////////////////////////////////////////
        startpos = part * para->getlimitOfNodesForVTK();
        endpos = startpos + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        for (uint i = 0; i < (uint)nodedatanames.size(); i++)
            nodedata[i].resize(sizeOfNodes);

        //////////////////////////////////////////////////////////////////////////
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {

                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsAreFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodedata[0][dn1] = (double)para->getParH(level)->pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[1][dn1] = (double)para->getParH(level)->rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodedata[2][dn1] = (double)para->getParH(level)->velocityX[pos] * (double)para->getVelocityRatio();
                nodedata[3][dn1] = (double)para->getParH(level)->velocityY[pos] * (double)para->getVelocityRatio();
                nodedata[4][dn1] = (double)para->getParH(level)->velocityZ[pos] * (double)para->getVelocityRatio();
                nodedata[5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];

                if(para->getDiffOn())
                    nodedata[firstConcNode][dn1] = (double)para->getParH(level)->concentration[pos];

                if(para->getIsBodyForce())
                {
                    nodedata[firstBodyForceNode    ][dn1] = (double)para->getParH(level)->forceX_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodedata[firstBodyForceNode + 1][dn1] = (double)para->getParH(level)->forceY_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodedata[firstBodyForceNode + 2][dn1] = (double)para->getParH(level)->forceZ_SP[pos] * (double)para->getScaledForceRatio(level);
                }

                if(para->getUseTurbulentViscosity())
                {
                    nodedata[firstNutNode][dn1] = (double)para->getParH(level)->turbViscosity[pos] * (double)para->getScaledViscosityRatio(level);
                }

                if (para->getCalcTurbulenceIntensity()) {
                    nodedata[firstTurbNode    ][dn1] = (double)para->getParH(level)->vxx[pos];
                    nodedata[firstTurbNode + 1][dn1] = (double)para->getParH(level)->vyy[pos];
                    nodedata[firstTurbNode + 2][dn1] = (double)para->getParH(level)->vzz[pos];
                    nodedata[firstTurbNode + 3][dn1] = (double)para->getParH(level)->vxy[pos];
                    nodedata[firstTurbNode + 4][dn1] = (double)para->getParH(level)->vxz[pos];
                    nodedata[firstTurbNode + 5][dn1] = (double)para->getParH(level)->vyz[pos];
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
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
            }
        }
        outFNames.push_back( WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata) );
    }
    return outFNames;
}

std::vector<std::string> FileWriter::writeUnstrucuredGridMedianLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
    std::vector< std::string > outFNames;

    std::vector< UbTupleFloat3 > nodes;
    std::vector< UbTupleUInt8 > cells;
    //std::vector< UbTupleUInt8 > cells2;
    std::vector< std::string > nodedatanames = getMedianNodeDataNames(para);
    int startIndex = para->getDiffOn()? 1 : 0;

    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    unsigned int dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsFluid;
    unsigned int startpos = 0;
    unsigned int endpos = 0;
    unsigned int sizeOfNodes = 0;
    std::vector< std::vector< double > > nodedata(nodedatanames.size());

    //printf("\n test for if... \n");
    for (unsigned int part = 0; part < fname.size(); part++)
    {
        //printf("\n test in if I... \n");
        //////////////////////////////////////////////////////////////////////////
        if (((part + 1) * para->getlimitOfNodesForVTK()) > (uint)para->getParH(level)->numberOfNodes)
        {
            sizeOfNodes = (uint)para->getParH(level)->numberOfNodes - (part * para->getlimitOfNodesForVTK());
        }
        else
        {
            sizeOfNodes = para->getlimitOfNodesForVTK();
        }
        //////////////////////////////////////////////////////////////////////////
        startpos = part * para->getlimitOfNodesForVTK();
        endpos = startpos + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        for (size_t i = 0; i < nodedatanames.size(); i++)
            nodedata[i].resize(sizeOfNodes);
        //////////////////////////////////////////////////////////////////////////
        //printf("\n test in if II... \n");
        for (unsigned int pos = startpos; pos < endpos; pos++)
        {
            if (para->getParH(level)->typeOfGridNode[pos] == GEO_FLUID)
            {
                //////////////////////////////////////////////////////////////////////////
                double x1 = para->getParH(level)->coordinateX[pos];
                double x2 = para->getParH(level)->coordinateY[pos];
                double x3 = para->getParH(level)->coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                number1 = pos;
                dn1 = pos - startpos;
                neighborsFluid = true;
                //////////////////////////////////////////////////////////////////////////
                nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                if(para->getDiffOn())
                    nodedata[0][dn1] = (double)para->getParH(level)->Conc_Med_Out[pos];
                nodedata[startIndex    ][dn1] = para->getParH(level)->press_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[startIndex + 1][dn1] = para->getParH(level)->rho_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
                nodedata[startIndex + 2][dn1] = para->getParH(level)->vx_SP_Med_Out[pos] * para->getVelocityRatio();
                nodedata[startIndex + 3][dn1] = para->getParH(level)->vy_SP_Med_Out[pos] * para->getVelocityRatio();
                nodedata[startIndex + 4][dn1] = para->getParH(level)->vz_SP_Med_Out[pos] * para->getVelocityRatio();
                nodedata[startIndex + 5][dn1] = (double)para->getParH(level)->typeOfGridNode[pos];
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
                if (number2 > endpos ||
                    number3 > endpos ||
                    number4 > endpos ||
                    number5 > endpos ||
                    number6 > endpos ||
                    number7 > endpos ||
                    number8 > endpos)  neighborsFluid = false;
                //////////////////////////////////////////////////////////////////////////
                dn2 = number2 - startpos;
                dn3 = number3 - startpos;
                dn4 = number4 - startpos;
                dn5 = number5 - startpos;
                dn6 = number6 - startpos;
                dn7 = number7 - startpos;
                dn8 = number8 - startpos;
                //////////////////////////////////////////////////////////////////////////
                if (isPeriodicCell(para, level, number2, number1, number3, number5))
                    continue;
                //////////////////////////////////////////////////////////////////////////
                if (neighborsFluid == true) cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
                //////////////////////////////////////////////////////////////////////////
            }
        }
        outFNames.push_back(WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata));
        //////////////////////////////////////////////////////////////////////////
    }
    return outFNames;
}