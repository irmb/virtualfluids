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
//! \file WriteMultiphaseQuantitiesCoProcessor.h
//! \ingroup CoProcessors
//! \author Konstantin Kutscher
//=======================================================================================
#include "WriteThixotropyQuantitiesCoProcessor.h"
#include "LBMKernel.h"
#include "BCProcessor.h"
#include "UbScheduler.h"
#include "DataSet3D.h"
#include "D3Q27System.h"
#include "BCArray3D.h"
#include <vector>
#include <string>
#include <algorithm> 
#include <numeric>
#include "basics/writer/WbWriterVtkXmlASCII.h"
#include "ThixotropyExpLBMKernel.h"
#include "Rheology.h"

using namespace std;

WriteThixotropyQuantitiesCoProcessor::WriteThixotropyQuantitiesCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
WriteThixotropyQuantitiesCoProcessor::WriteThixotropyQuantitiesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string& path, WbWriter* const writer, SPtr<LBMUnitConverter> conv, std::shared_ptr<vf::mpi::Communicator> comm) : CoProcessor(grid, s), path(path), writer(writer),	conv(conv),	comm(comm)
{
	gridRank = comm->getProcessID();
	minInitLevel = this->grid->getCoarsestInitializedLevel();
	maxInitLevel = this->grid->getFinestInitializedLevel();

	blockVector.resize(maxInitLevel + 1);

	for (int level = minInitLevel; level <= maxInitLevel; level++)
	{
		grid->getBlocks(level, gridRank, true, blockVector[level]);
	}
}
//////////////////////////////////////////////////////////////////////////
void WriteThixotropyQuantitiesCoProcessor::init()
{

}
//////////////////////////////////////////////////////////////////////////
void WriteThixotropyQuantitiesCoProcessor::process(double step)
{
	if (scheduler->isDue(step))
		collectData(step);

	UBLOG(logDEBUG3, "WriteThixotropyQuantitiesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void WriteThixotropyQuantitiesCoProcessor::collectData(double step)
{
	int istep = static_cast<int>(step);
	//ConcentrationSum = 0;
	for (int level = minInitLevel; level <= maxInitLevel; level++)
	{
		for(SPtr<Block3D> block : blockVector[level])
		{
			if (block)
			{
				addDataMQ(block);
			}
		}
	}
	//cout <<"at step = "<< step << ", the sum of Concentration in the domain = "<<ConcentrationSum << endl;
	string pfilePath, partPath, subfolder, cfilePath;

	subfolder = "thyxotropy" + UbSystem::toString(istep);
	pfilePath = path + "/thyxotropy/" + subfolder;
	cfilePath = path + "/thyxotropy/thyxotropy_collection";
	partPath = pfilePath + "/thyxotropy" + UbSystem::toString(gridRank) + "_" + UbSystem::toString(istep);


	string partName = writer->writeOctsWithNodeData(partPath, nodes, cells, datanames, data);

	size_t found = partName.find_last_of("/");
	string piece = partName.substr(found + 1);
	piece = subfolder + "/" + piece;

	vector<string> cellDataNames;
	vector<string> pieces = comm->gather(piece);
	if (comm->getProcessID() == comm->getRoot())
	{
		string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath, pieces, datanames, cellDataNames);
		found = pname.find_last_of("/");
		piece = pname.substr(found + 1);

		vector<string> filenames;
		filenames.push_back(piece);
		if (step == CoProcessor::scheduler->getMinBegin())
		{
			WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath, filenames, istep, false);
		}
		else
		{
			WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath, filenames, istep, false);
		}
		UBLOG(logINFO, "WriteThixotropyQuantitiesCoProcessor step: " << istep);
	}

	clearData();
}
//////////////////////////////////////////////////////////////////////////
void WriteThixotropyQuantitiesCoProcessor::clearData()
{
	nodes.clear();
	cells.clear();
	datanames.clear();
	data.clear();
}
//////////////////////////////////////////////////////////////////////////
void WriteThixotropyQuantitiesCoProcessor::addDataMQ(SPtr<Block3D> block)
{
	UbTupleDouble3 org = grid->getBlockWorldCoordinates(block);;
	UbTupleDouble3 nodeOffset = grid->getNodeOffset(block);
	double         dx = grid->getDeltaX(block);

	//double level = (double)block->getLevel();
	//double blockID = (double)block->getGlobalID();

	//Diese Daten werden geschrieben:
	datanames.resize(0);
	datanames.push_back("viscosity");
	//datanames.push_back("lambda");
	//datanames.push_back("ShearRate");
	datanames.push_back("omega");
	//datanames.push_back("Fluxx");
	//datanames.push_back("Fluxy");
	//datanames.push_back("Fluxz");
	//datanames.push_back("delC");
	//datanames.push_back("diff");
	//datanames.push_back("AccumulatedConcentration");
	//datanames.push_back("m100");
	//datanames.push_back("Level");
	//datanames.push_back("BlockID");



	data.resize(datanames.size());

   SPtr<ILBMKernel> kernel = block->getKernel();
   SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();          
   SPtr<DistributionArray3D> distributionsF = kernel->getDataSet()->getFdistributions(); 
	//SPtr<DistributionArray3D> distributionsH = kernel->getDataSet()->getHdistributions();
	//LBMReal collFactorF = staticPointerCast<ThixotropyExpLBMKernel>(kernel)->getCollisionFactorF();
	LBMReal collFactor = kernel->getCollisionFactor();
	LBMReal f[D3Q27System::ENDF + 1];
	//LBMReal h[D3Q27System::ENDF + 1];
	//LBMReal viscosity=0; // lambda, gammaDot;
	
	//knotennummerierung faengt immer bei 0 an!
	int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

	int minX1 = 0;
	int minX2 = 0;
	int minX3 = 0;

	int maxX1 = (int)(distributionsF->getNX1());
	int maxX2 = (int)(distributionsF->getNX2());
	int maxX3 = (int)(distributionsF->getNX3());

	//int minX1 = 1;
	//int minX2 = 1;
	//int minX3 = 1;

	//int maxX1 = (int)(distributions->getNX1());
	//int maxX2 = (int)(distributions->getNX2());
	//int maxX3 = (int)(distributions->getNX3());

	//nummern vergeben und node vector erstellen + daten sammeln
	CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
	maxX1 -= 2;
	maxX2 -= 2;
	maxX3 -= 2;
	//D3Q27BoundaryConditionPtr bcPtr;
	int nr = (int)nodes.size();

	for (int ix3 = minX3; ix3 <= maxX3; ix3++)
	{
		for (int ix2 = minX2; ix2 <= maxX2; ix2++)
		{
			for (int ix1 = minX1; ix1 <= maxX1; ix1++)
			{
				if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3))
				{
					int index = 0;
					nodeNumbers(ix1, ix2, ix3) = nr++;
					nodes.push_back(makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
						float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
						float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)));

					//distributionsH->getDistribution(h, ix1, ix2, ix3);

					//lambda = D3Q27System::getDensity(h);

					//if (UbMath::isNaN(lambda) || UbMath::isInfinity(lambda) )
					//	UB_THROW(UbException(UB_EXARGS, "lambda is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" + block->toString() +
					//		", node=" + UbSystem::toString(ix1) + "," + UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));

					//distributionsF->getDistribution(f, ix1, ix2, ix3);
					//LBMReal rho = D3Q27System::getDensity(f);
					
					//LBMReal gammaDot = D3Q27System::getShearRate(f, collFactor);

					//LBMReal collFactorF = collFactor - 1e-6 / (gammaDot + one * 1e-9);
					//collFactorF = (collFactorF < 0.5) ? 0.5 : collFactorF;

					//LBMReal collFactorF = Rheology::getBinghamCollFactor(collFactor, gammaDot, rho);

					//data[index++].push_back(lambda);
					//data[index++].push_back(gammaDot);
					//data[index++].push_back(collFactorF);

					distributionsF->getDistribution(f, ix1, ix2, ix3);
					LBMReal rho = D3Q27System::getDensity(f);
					LBMReal shearRate = D3Q27System::getShearRate(f, collFactor);
					//LBMReal omega = Rheology::getHerschelBulkleyCollFactor(collFactor, shearRate, rho);
					//LBMReal omega = Rheology::getPowellEyringCollFactor(collFactor, shearRate, rho);
					LBMReal omega = Rheology::getBinghamCollFactor(collFactor, shearRate, rho);
					LBMReal viscosity = (omega == 0) ? 0 : UbMath::c1o3 * (UbMath::c1/omega-UbMath::c1o2);

					
					data[index++].push_back(viscosity);
					data[index++].push_back(omega);
				}
			}
		}
	}
	maxX1 -= 1;
	maxX2 -= 1;
	maxX3 -= 1;
	//cell vector erstellen
	for (int ix3 = minX3; ix3 <= maxX3; ix3++)
	{
		for (int ix2 = minX2; ix2 <= maxX2; ix2++)
		{
			for (int ix1 = minX1; ix1 <= maxX1; ix1++)
			{
				if ((SWB = nodeNumbers(ix1, ix2, ix3)) >= 0
					&& (SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0
					&& (NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0
					&& (NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0
					&& (SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0
					&& (SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0
					&& (NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0
					&& (NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0)
				{
					cells.push_back(makeUbTuple((unsigned int)SWB, (unsigned int)SEB, (unsigned int)NEB, (unsigned int)NWB, (unsigned int)SWT, (unsigned int)SET, (unsigned int)NET, (unsigned int)NWT));
				}
			}
		}
	}
}
