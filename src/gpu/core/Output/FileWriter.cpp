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

#include <basics/StringUtilities/StringUtil.h>
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"
#include "WriterUtilities.h"

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"

std::string makeCollectionFileNameEnding(int ID, int timestep)
{
    return "_ID_" + StringUtil::toString<int>(ID) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";
}

std::string makePartFileName(const std::string &prefix, uint level, int ID, int part, int timestep)
{
    return prefix + "_bin" + WriterUtilities::makePartFileNameEnding(level, ID, part, timestep);
}

std::string makeMedianPartFileName(const std::string &prefix, uint level, int ID, int part, int timestep)
{
    return prefix + "_bin_median" + WriterUtilities::makePartFileNameEnding(level, ID, part, timestep);
}


std::string makeCollectionFileName(const std::string &prefix, int ID, int timestep)
{
    return prefix + "_bin" + makeCollectionFileNameEnding(ID, timestep);
}

std::string makeMedianCollectionFileName(const std::string &prefix, int ID, int timestep)
{
    return prefix + "_bin_median" + makeCollectionFileNameEnding(ID, timestep);
}

std::string makePvdCollectionFileName(const std::string &prefix, int mpiProcessID)
{
    return prefix + "_bin" + "_ID_" + StringUtil::toString<int>(mpiProcessID);
}


void FileWriter::writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaMemoryManager)
{
    unsigned int timestep = para->getTimestepInit();
    for (int level = para->getCoarse(); level <= para->getFine(); level++) {
        cudaMemoryManager->cudaCopyPrint(level);
        writeTimestep(para, timestep, level);
    }

    this->fileNamesForCollectionFileTimeSeries[timestep] = this->fileNamesForCollectionFile;
    this->writeCollectionFile(para, timestep);

    if( para->getCalcMedian() )
        this->writeCollectionFileMedian(para, timestep);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    for (int level = para->getCoarse(); level <= para->getFine(); level++)
        writeTimestep(para, timestep, level);

    this->fileNamesForCollectionFileTimeSeries[timestep] = this->fileNamesForCollectionFile;
    if (timestep == para->getTimestepEnd())
        this->writePvdCollectionFileForTimeSeries(*para);

    this->writeCollectionFile(para, timestep);

    if( para->getCalcMedian() )
        this->writeCollectionFileMedian(para, timestep);
}

void FileWriter::writeTimestep(std::shared_ptr<Parameter> para, unsigned int timestep, int level)
{
    const unsigned int numberOfParts = WriterUtilities::calculateNumberOfParts(para.get(), level);
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

std::string FileWriter::writeCollectionFile(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    auto nodeDataNames = this->getNodeDataNames(para);
    std::vector<std::string> cellDataNames;
    std::string pFileName= WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFile, nodeDataNames, cellDataNames);
    this->fileNamesForCollectionFile.clear();
    return pFileName;
}

std::string FileWriter::writeCollectionFileMedian(std::shared_ptr<Parameter> para, unsigned int timestep)
{
    std::string filename = makeMedianCollectionFileName(para->getFName(), para->getMyProcessID(), timestep);
    std::vector<std::string> nodeDataNames = getMedianNodeDataNames(para);
    std::vector<std::string> cellDataNames;
    std::string pFileName =  WbWriterVtkXmlBinary::getInstance()->writeParallelFile(filename, this->fileNamesForCollectionFileMedian, nodeDataNames, cellDataNames);
    this->fileNamesForCollectionFileMedian.clear();
    return pFileName;
}

std::string FileWriter::writePvdCollectionFileForTimeSeries(const Parameter &para)
{
    const std::string filename = makePvdCollectionFileName(para.getFName(), para.getMyProcessID());
    return WbWriterVtkXmlBinary::getInstance()->writeCollectionForTimeSeries(filename, this->fileNamesForCollectionFileTimeSeries, false);
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

    std::array<uint, 8> indicesOfOct;
    std::array<uint, 8> relativePosInPart;
    uint relPosInPart;
    bool allNodesValid;
    unsigned int startPosition;
    unsigned int endPosition;
    unsigned int sizeOfNodes;
    std::vector< std::vector< double > > nodeData(nodeDataNames.size());

    for (unsigned int part = 0; part < fname.size(); part++)
    {
        sizeOfNodes = WriterUtilities::calculateNumberOfNodesInPart(para.get(), level, part);

        //////////////////////////////////////////////////////////////////////////
        startPosition = part * para->getLimitOfNodesForVTK();
        endPosition = startPosition + sizeOfNodes;
        //////////////////////////////////////////////////////////////////////////
        cells.clear();
        nodes.resize(sizeOfNodes);
        for (uint i = 0; i < (uint)nodeDataNames.size(); i++)
            nodeData[i].resize(sizeOfNodes);

        //////////////////////////////////////////////////////////////////////////
        for (unsigned int pos = startPosition; pos < endPosition; pos++) {
            const LBMSimulationParameter& parH = para->getParHostAsReference(level);
            if (parH.typeOfGridNode[pos] == GEO_FLUID) {
                //////////////////////////////////////////////////////////////////////////
                double x1 = parH.coordinateX[pos];
                double x2 = parH.coordinateY[pos];
                double x3 = parH.coordinateZ[pos];
                //////////////////////////////////////////////////////////////////////////
                relPosInPart = pos - startPosition;
                //////////////////////////////////////////////////////////////////////////
                nodes[relPosInPart] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
                nodeData[0][relPosInPart] = (double)parH.pressure[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodeData[1][relPosInPart] = (double)parH.rho[pos] / (double)3.0 * (double)para->getDensityRatio() * (double)para->getVelocityRatio() * (double)para->getVelocityRatio();
                nodeData[2][relPosInPart] = (double)parH.velocityX[pos] * (double)para->getVelocityRatio();
                nodeData[3][relPosInPart] = (double)parH.velocityY[pos] * (double)para->getVelocityRatio();
                nodeData[4][relPosInPart] = (double)parH.velocityZ[pos] * (double)para->getVelocityRatio();
                nodeData[5][relPosInPart] = (double)parH.typeOfGridNode[pos];

                if(para->getDiffOn())
                    nodeData[firstConcNode][relPosInPart] = (double)parH.concentration[pos];

                if(para->getIsBodyForce())
                {
                    nodeData[firstBodyForceNode    ][relPosInPart] = (double)parH.forceX_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodeData[firstBodyForceNode + 1][relPosInPart] = (double)parH.forceY_SP[pos] * (double)para->getScaledForceRatio(level);
                    nodeData[firstBodyForceNode + 2][relPosInPart] = (double)parH.forceZ_SP[pos] * (double)para->getScaledForceRatio(level);
                }

                if(para->getUseTurbulentViscosity())
                {
                    nodeData[firstNutNode][relPosInPart] = (double)parH.turbViscosity[pos] * (double)para->getScaledViscosityRatio(level);
                }

                if (para->getCalcTurbulenceIntensity()) {
                    nodeData[firstTurbulenceNode    ][relPosInPart] = (double)parH.vxx[pos];
                    nodeData[firstTurbulenceNode + 1][relPosInPart] = (double)parH.vyy[pos];
                    nodeData[firstTurbulenceNode + 2][relPosInPart] = (double)parH.vzz[pos];
                    nodeData[firstTurbulenceNode + 3][relPosInPart] = (double)parH.vxy[pos];
                    nodeData[firstTurbulenceNode + 4][relPosInPart] = (double)parH.vxz[pos];
                    nodeData[firstTurbulenceNode + 5][relPosInPart] = (double)parH.vyz[pos];
                }

                //////////////////////////////////////////////////////////////////////////

                WriterUtilities::getIndicesOfAllNodesInOct(indicesOfOct, pos, para->getParHostAsReference(level));

                if (WriterUtilities::isPeriodicCell(parH, indicesOfOct[0], indicesOfOct[6])) {
                    continue;
                }

                //////////////////////////////////////////////////////////////////////////
                allNodesValid = WriterUtilities::areAllNodesInOctValidForWriting(indicesOfOct, parH, endPosition);

                //////////////////////////////////////////////////////////////////////////
                if (allNodesValid) {
                    WriterUtilities::calculateRelativeNodeIndexInPart(relativePosInPart, indicesOfOct, startPosition);
                    cells.push_back(makeUbTupleFromArray(relativePosInPart));
                }
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
        sizeOfNodes = WriterUtilities::calculateNumberOfNodesInPart(para.get(), level, part);
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
                if (WriterUtilities::isPeriodicCell(para->getParHostAsReference(level), number1, number7))
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