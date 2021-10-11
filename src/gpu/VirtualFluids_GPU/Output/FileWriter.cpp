//  _    ___      __              __________      _     __        ______________   __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____   /  ___/ __  / /  / /
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/  / /___/ /_/ / /  / /
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )  / /_) / ____/ /__/ / 
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/   \____/_/    \_____/
//
//////////////////////////////////////////////////////////////////////////
#include "FileWriter.h"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <cmath>

#include <Core/StringUtilities/StringUtil.h>

#include "Parameter/Parameter.h"
#include "GPU/CudaMemoryManager.h"

#include "LBM/LB.h"
#include "LBM/D3Q27.h"

#include <basics/writer/WbWriterVtkXmlBinary.h>

void FileWriter::writeInit(std::shared_ptr<Parameter> para, std::shared_ptr<CudaMemoryManager> cudaManager)
{
    unsigned int timestep = para->getTInit();
	for (int level = para->getCoarse(); level <= para->getFine(); level++) {
		cudaManager->cudaCopyPrint(level);
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
    const unsigned int numberOfParts = para->getParH(level)->size_Mat_SP / para->getlimitOfNodesForVTK() + 1;
    std::vector<std::string> fname;
	std::vector<std::string> fnameMed;
	for (unsigned int i = 1; i <= numberOfParts; i++)
	{
		fname.push_back(para->getFName() + "_bin_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(i) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk");
		fnameMed.push_back(para->getFName() + "_bin_median_lev_" + StringUtil::toString<int>(level) + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_Part_" + StringUtil::toString<int>(i) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk");

        this->fileNamesForCollectionFile.push_back( fname.back() );
        this->fileNamesForCollectionFileMedian.push_back( fnameMed.back() );
	}

	if (para->getDiffOn() == true)
		writeUnstrucuredGridLTConc(para, level, fname);
	else
		writeUnstrucuredGridLT(para, level, fname);

	if (para->getCalcMedian())
	{
		if (para->getDiffOn() == true)
			writeUnstrucuredGridMedianLTConc(para, level, fnameMed);
		else
			writeUnstrucuredGridMedianLT(para, level, fnameMed);
	}
}

bool FileWriter::isPeriodicCell(std::shared_ptr<Parameter> para, int level, unsigned int number2, unsigned int number1, unsigned int number3, unsigned int number5)
{
	return (para->getParH(level)->coordX_SP[number2] < para->getParH(level)->coordX_SP[number1]) ||
		   (para->getParH(level)->coordY_SP[number3] < para->getParH(level)->coordY_SP[number1]) ||
		   (para->getParH(level)->coordZ_SP[number5] < para->getParH(level)->coordZ_SP[number1]);
}

void VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFile(std::shared_ptr<Parameter> para, unsigned int timestep)
{

    std::string filename = para->getFName() + "_bin_ID_" + StringUtil::toString<int>(para->getMyID()) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";

    std::ofstream file;

    file.open( filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PPointData>" << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"press\" /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"rho\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx1\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx2\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx3\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"geo\"   /> " << std::endl;
    if( para->getDiffOn() ) file << "       <DataArray type=\"Float32\" Name=\"conc\"  /> " << std::endl;
    file << "    </PPointData>" << std::endl;

    file << "    <PPoints>" << std::endl;
    file << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << std::endl;
    file << "    </PPoints>" << std::endl;

    for( auto& fname : this->fileNamesForCollectionFile )
    {
        const auto filenameWithoutPath=fname.substr( fname.find_last_of('/') + 1 );
        file << "    <Piece Source=\"" << filenameWithoutPath << ".bin.vtu\"/>" << std::endl;
    }

    file << "  </PUnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file.close();

    this->fileNamesForCollectionFile.clear();
}

void VIRTUALFLUIDS_GPU_EXPORT FileWriter::writeCollectionFileMedian(std::shared_ptr<Parameter> para, unsigned int timestep)
{

    std::string filename = para->getFName() + "_bin_median_ID_" + StringUtil::toString<int>(para->getMyID()) + "_t_" + StringUtil::toString<int>(timestep) + ".vtk";

    std::ofstream file;

    file.open( filename + ".pvtu" );

    //////////////////////////////////////////////////////////////////////////
    
    file << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
    file << "  <PUnstructuredGrid GhostLevel=\"1\">" << std::endl;

    file << "    <PPointData>" << std::endl;
    if( para->getDiffOn() ) file << "       <DataArray type=\"Float32\" Name=\"concMed\"  /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"pressMed\" /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"rhoMed\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx1Med\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx2Med\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"vx3Med\"   /> " << std::endl;
    file << "       <DataArray type=\"Float32\" Name=\"geo\"   /> " << std::endl;
    file << "    </PPointData>" << std::endl;

    file << "    <PPoints>" << std::endl;
    file << "      <PDataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\"/>" << std::endl;
    file << "    </PPoints>" << std::endl;

    for( auto& fname : this->fileNamesForCollectionFileMedian )
    {
        const auto filenameWithoutPath=fname.substr( fname.find_last_of('/') + 1 );
        file << "    <Piece Source=\"" << filenameWithoutPath << ".bin.vtu\"/>" << std::endl;
    }

    file << "  </PUnstructuredGrid>" << std::endl;
    file << "</VTKFile>" << std::endl;

    //////////////////////////////////////////////////////////////////////////

    file.close();

    this->fileNamesForCollectionFileMedian.clear();
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
    //nodedatanames.push_back("sendNodes");
    //nodedatanames.push_back("sparseIndex");
    unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
    uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
    bool neighborsAreFluid;
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
        for (uint i = 0; i < (uint)nodedatanames.size(); i++)
			nodedata[i].resize(sizeOfNodes);

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

				//int sendNode = 0; // 0 - not a sendNode; 1 - sendNode; 2 - sendNode in communication after fine to coarse
    //            testForSendNodeY(para, level, pos, sendNode); // slow and should not be done multiple times --> use for debugging only!
				//nodedata[6][dn1] = (double) sendNode;
    //            nodedata[7][dn1] = (double) pos;
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
				if (isPeriodicCell(para, level, number2, number1, number3, number5))
					continue;
				//////////////////////////////////////////////////////////////////////////
                if (neighborsAreFluid)
                    cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
            }
        }
        WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
    }
}

void FileWriter::testForSendNodeY(std::shared_ptr<Parameter> &para, int level, unsigned int pos, int &sendNode)
{
    testForCommunicationNode(para, level, pos, sendNode, &para->getParH(level)->sendProcessNeighborY,
                        &para->getParH(level)->sendProcessNeighborsAfterFtoCY);
}

void FileWriter::testForCommunicationNode(std::shared_ptr<Parameter> &para, int level, unsigned int pos, int &sendNode,
                                     std::vector<PN27> *sendOrRecvProcessNeighbor,
                                     std::vector<PN27> *sendOrRecvProcessNeighborAfterFtoC)
{
    for (uint direction = 0; direction < (uint)sendOrRecvProcessNeighbor->size(); direction++) {
        for (int i = 0; i < (*sendOrRecvProcessNeighbor)[direction].numberOfNodes; i++) {
            if (pos == (uint)(*sendOrRecvProcessNeighbor)[direction].index[i]) {
                sendNode = 1;
				//std::cout << "send/recv Node found: " << pos << std::endl;
                if (level < para->getMaxLevel() && i < (*sendOrRecvProcessNeighborAfterFtoC)[direction].numberOfNodes) {
                    //std::cout << "Communication send/recv Node found: " << pos << std::endl;
					sendNode = 2;
                }
                return;
            }
        }
    }
}

void FileWriter::writeUnstrucuredGridLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
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
	nodedatanames.push_back("conc");
	unsigned int number1, number2, number3, number4, number5, number6, number7, number8;
	uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
	bool neighborsAreFluid;
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
		nodedata[6].resize(sizeOfNodes);
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
				nodedata[6][dn1] = (double)para->getParH(level)->Conc[pos];
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
				if (isPeriodicCell(para, level, number2, number1, number3, number5))
					continue;
				//////////////////////////////////////////////////////////////////////////
				if (neighborsAreFluid)
					cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
			}
		}
		WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
	}
}

void FileWriter::writeUnstrucuredGridMedianLT(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
	std::vector< UbTupleFloat3 > nodes;
	std::vector< UbTupleUInt8 > cells;
	//std::vector< UbTupleUInt8 > cells2;
	std::vector< std::string > nodedatanames;
	nodedatanames.push_back("pressMed");
	nodedatanames.push_back("rhoMed");
	nodedatanames.push_back("vx1Med");
	nodedatanames.push_back("vx2Med");
	nodedatanames.push_back("vx3Med");
	nodedatanames.push_back("geo");
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
		if (((part + 1)*para->getlimitOfNodesForVTK()) > para->getParH(level)->size_Mat_SP)
		{
			sizeOfNodes = para->getParH(level)->size_Mat_SP - (part * para->getlimitOfNodesForVTK());
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
		nodedata[0].resize(sizeOfNodes);
		nodedata[1].resize(sizeOfNodes);
		nodedata[2].resize(sizeOfNodes);
		nodedata[3].resize(sizeOfNodes);
		nodedata[4].resize(sizeOfNodes);
		nodedata[5].resize(sizeOfNodes);
		//////////////////////////////////////////////////////////////////////////
		//printf("\n test in if II... \n");
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
				neighborsFluid = true;
				//////////////////////////////////////////////////////////////////////////
				nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
				nodedata[0][dn1] = para->getParH(level)->press_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
				nodedata[1][dn1] = para->getParH(level)->rho_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
				nodedata[2][dn1] = para->getParH(level)->vx_SP_Med_Out[pos] * para->getVelocityRatio();
				nodedata[3][dn1] = para->getParH(level)->vy_SP_Med_Out[pos] * para->getVelocityRatio();
				nodedata[4][dn1] = para->getParH(level)->vz_SP_Med_Out[pos] * para->getVelocityRatio();
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
					para->getParH(level)->geoSP[number8] != GEO_FLUID)  neighborsFluid = false;
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
		WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
		//////////////////////////////////////////////////////////////////////////
	}
}

void FileWriter::writeUnstrucuredGridMedianLTConc(std::shared_ptr<Parameter> para, int level, std::vector<std::string >& fname)
{
	std::vector< UbTupleFloat3 > nodes;
	std::vector< UbTupleUInt8 > cells;
	std::vector< std::string > nodedatanames;
	nodedatanames.push_back("concMed");
	nodedatanames.push_back("pressMed");
	nodedatanames.push_back("rhoMed");
	nodedatanames.push_back("vx1Med");
	nodedatanames.push_back("vx2Med");
	nodedatanames.push_back("vx3Med");
	nodedatanames.push_back("geo");
	uint number1, number2, number3, number4, number5, number6, number7, number8;
	uint dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8;
	bool neighborsFluid;
	uint startpos = 0;
	uint endpos = 0;
	uint sizeOfNodes = 0;
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
		nodedata[6].resize(sizeOfNodes);
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
				neighborsFluid = true;
				//////////////////////////////////////////////////////////////////////////
				nodes[dn1] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
				nodedata[0][dn1] = (double)para->getParH(level)->Conc_Med_Out[pos];
				nodedata[1][dn1] = (double)para->getParH(level)->press_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
				nodedata[2][dn1] = (double)para->getParH(level)->rho_SP_Med_Out[pos] / 3.0f * para->getDensityRatio() * para->getVelocityRatio() * para->getVelocityRatio();
				nodedata[3][dn1] = (double)para->getParH(level)->vx_SP_Med_Out[pos] * para->getVelocityRatio();
				nodedata[4][dn1] = (double)para->getParH(level)->vy_SP_Med_Out[pos] * para->getVelocityRatio();
				nodedata[5][dn1] = (double)para->getParH(level)->vz_SP_Med_Out[pos] * para->getVelocityRatio();
				nodedata[6][dn1] = (double)para->getParH(level)->geoSP[pos];
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
					para->getParH(level)->geoSP[number8] != GEO_FLUID)  neighborsFluid = false;
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
				if (neighborsFluid) 
					cells.push_back(makeUbTuple(dn1, dn2, dn3, dn4, dn5, dn6, dn7, dn8));
				//////////////////////////////////////////////////////////////////////////
			}
		}
		WbWriterVtkXmlBinary::getInstance()->writeOctsWithNodeData(fname[part], nodes, cells, nodedatanames, nodedata);
		//////////////////////////////////////////////////////////////////////////
	}
}
//////////////////////////////////////////////////////////////////////////









