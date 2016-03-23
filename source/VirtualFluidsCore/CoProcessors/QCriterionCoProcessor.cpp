#include "QCriterionCoProcessor.h"
#include "LBMKernel3D.h"
#include "LBMKernelETD3Q27.h"
#include "D3Q27ETBCProcessor.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>
#include "basics/writer/WbWriterVtkXmlASCII.h"


using namespace std;
using namespace D3Q27System;

QCriterionCoProcessor::QCriterionCoProcessor(Grid3DPtr grid, const std::string& path, 
	WbWriter* const writer,
	UbSchedulerPtr s, CommunicatorPtr comm)
	: CoProcessor(grid, s),
	path(path),
	comm(comm),
	writer(writer)
{
	init();
}
//////////////////////////////////////////////////////////////////////////
void QCriterionCoProcessor::init()
{
	gridRank  = comm->getProcessID(); 
	minInitLevel = this->grid->getCoarsestInitializedLevel();
	maxInitLevel = this->grid->getFinestInitializedLevel();

	blockVector.resize(maxInitLevel+1);

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		grid->getBlocks(level, gridRank, true, blockVector[level]); //grid: private variable in Postprocessor. Initialized by filling with blocks
	}
}
//////////////////////////////////////////////////////////////////////////
void QCriterionCoProcessor::process(double step)
{
	if(scheduler->isDue(step) )
		collectData(step);

	UBLOG(logDEBUG3, "QCriterionPostprocessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void QCriterionCoProcessor::collectData(double step)
{
	int istep = static_cast<int>(step);

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		BOOST_FOREACH(Block3DPtr block, blockVector[level])
		{
			if (block)
			{
				addData(block);

			}
		}
	}

	string partName = writer->writeOctsWithNodeData(path+ UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep),nodes,cells,datanames,data);
	size_t found=partName.find_last_of("//");
	string piece = partName.substr(found+1);

	vector<string> cellDataNames;

	//distributed writing as in MacroscopicValuesPostprocessor.cpp
	vector<string> pieces = comm->gather(piece); //comm: MPI-Wrapper
	if (comm->getProcessID() == comm->getRoot())
	{
		string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(path+"_"+UbSystem::toString(istep),pieces,datanames,cellDataNames);

		vector<string> filenames;
		filenames.push_back(pname);
		if (step == CoProcessor::scheduler->getMinBegin()) //first time in timeseries
		{
			WbWriterVtkXmlASCII::getInstance()->writeCollection(path+"_collection",filenames,istep,false);
		} 
		else
		{
			WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(path+"_collection",filenames,istep,false);
		}
		UBLOG(logINFO,"QCriterionPostprocessor step: " << istep);
	}

	clearData();


}
//////////////////////////////////////////////////////////////////////////
void QCriterionCoProcessor::clearData()
{
	nodes.clear();
	cells.clear();
	datanames.clear();
	data.clear();
}
//////////////////////////////////////////////////////////////////////////
void QCriterionCoProcessor::addData(const Block3DPtr block)
{
	UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
	UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
	UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
	double         dx           = grid->getDeltaX(block);

	//Diese Daten werden geschrieben:
	datanames.resize(0);
	datanames.push_back("q");
	datanames.push_back("scaleFactor");
	data.resize(datanames.size());


	LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
	BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();          
	DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();     

	int SWB,SEB,NEB,NWB,SWT,SET,NET,NWT;

	int minX1 = 0;
	int minX2 = 0;
	int minX3 = 0;

	int maxX1 = (int)(distributions->getNX1());
	int maxX2 = (int)(distributions->getNX2());
	int maxX3 = (int)(distributions->getNX3());

	int currentLevel = block->getLevel();
	//nummern vergeben und node vector erstellen + daten sammeln
	CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3,-1);
	maxX1 -= 2; //-2 wegen ghost layer: 
	maxX2 -= 2; //0-maxXi-1 ist arraygroesse. 
	maxX3 -= 2; //ueberlapp 1 in +,- Richtung. zum schreiben werden statt feldern von 1 bis (max-2) felder von 0 bis max-3 verwendet! 


	int nr = (int)nodes.size();

	for(int ix3=minX3; ix3<=maxX3; ix3++)
	{
		for(int ix2=minX2; ix2<=maxX2; ix2++)
		{
			for(int ix1=minX1; ix1<=maxX1; ix1++)
			{
				if(!bcArray.isUndefined(ix1,ix2,ix3) && !bcArray.isSolid(ix1,ix2,ix3))
				{
					//nodeNumbers-vektor wird mit koordinaten befuellt
					int index = 0;
					nodeNumbers(ix1,ix2,ix3) = nr++;
					nodes.push_back( makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
						float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
						float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)) );			

					/////////////////////////////
					// Geschwindigkeitsvektoren 
					LBMReal vE[3];
					LBMReal vW[3];
					LBMReal vN[3];
					LBMReal vS[3];
					LBMReal vT[3];
					LBMReal vB[3];
					//hole geschwindigkeiten an nachbarknoten
					getNeighborVelocities(1,0,0, ix1,  ix2,  ix3, block, vE, vW);
					getNeighborVelocities(0,1,0, ix1,  ix2,  ix3, block, vN, vS);
					getNeighborVelocities(0,0,1, ix1,  ix2,  ix3, block, vT, vB);
					//////////////////////////////////
					//derivatives
					LBMReal duxdy=(vN[xdir]-vS[xdir])*0.5;
					LBMReal duydx=(vE[ydir]-vW[ydir])*0.5;
					LBMReal duxdz=(vT[xdir]-vB[xdir])*0.5;
					LBMReal duzdx=(vE[zdir]-vW[zdir])*0.5;
					LBMReal duydz=(vT[ydir]-vB[ydir])*0.5;
					LBMReal duzdy=(vN[zdir]-vS[zdir])*0.5;

					LBMReal duxdx=(vE[xdir]-vW[xdir])*0.5;
					LBMReal duydy=(vN[ydir]-vS[ydir])*0.5;
					LBMReal duzdz=(vT[zdir]-vB[zdir])*0.5;					

					LBMReal scaleFactor=(double)(1<<(currentLevel-minInitLevel));//pow(2.0,(double)(currentLevel-minInitLevel));//finer grid -> current level higher. coarsest grid: currentLevel=minInitLevel=0
					// Q=-0.5*(S_ij S_ij - Omega_ij Omega_ij) => regions where vorticity is larger than strain rate
					LBMReal q=-(duxdy*duydx+duxdz*duzdx+duydz*duzdy+duxdx*duxdx+duydy*duydy+duzdz*duzdz)*scaleFactor;

					data[index++].push_back( q );
					data[index++].push_back( scaleFactor );

				}
			}
		}
	}
	maxX1 -= 1;
	maxX2 -= 1;
	maxX3 -= 1;
	//cell vector erstellen
	for(int ix3=minX3; ix3<=maxX3; ix3++)
	{
		for(int ix2=minX2; ix2<=maxX2; ix2++)
		{
			for(int ix1=minX1; ix1<=maxX1; ix1++)
			{
				if(   (SWB=nodeNumbers( ix1  , ix2,   ix3   )) >= 0
					&& (SEB=nodeNumbers( ix1+1, ix2,   ix3   )) >= 0
					&& (NEB=nodeNumbers( ix1+1, ix2+1, ix3   )) >= 0
					&& (NWB=nodeNumbers( ix1  , ix2+1, ix3   )) >= 0 
					&& (SWT=nodeNumbers( ix1  , ix2,   ix3+1 )) >= 0
					&& (SET=nodeNumbers( ix1+1, ix2,   ix3+1 )) >= 0
					&& (NET=nodeNumbers( ix1+1, ix2+1, ix3+1 )) >= 0
					&& (NWT=nodeNumbers( ix1  , ix2+1, ix3+1 )) >= 0                )
				{
					// for valid points: neighbors are added to cells-vector
					cells.push_back( makeUbTuple(SWB,SEB,NEB,NWB,SWT,SET,NET,NWT) );
				}
			}
		}
	}

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QCriterionCoProcessor::getNeighborVelocities(int offx, int offy, int offz, int ix1, int ix2, int ix3, const Block3DPtr block, LBMReal* vE, LBMReal* vW)
{
	LBMKernelETD3Q27Ptr kernel = boost::dynamic_pointer_cast<LBMKernelETD3Q27>(block->getKernel());
	BCArray3D<D3Q27BoundaryCondition>& bcArray = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernel->getBCProcessor())->getBCArray();          
	DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();   

	int minX1 = 0;
	int minX2 = 0;
	int minX3 = 0;

	int maxX1 = (int)(distributions->getNX1());
	int maxX2 = (int)(distributions->getNX2());
	int maxX3 = (int)(distributions->getNX3());
	if (maxX1<3) throw UbException(UB_EXARGS,"QCriterionPostprocessor: NX1 too small for FD stencils!");
	maxX1 -= 2;
	maxX2 -= 2;
	maxX3 -= 2;
	bool checkInterpolation=true;
	bool neighNodeIsBC=false;
	D3Q27BoundaryConditionPtr bcPtr;

	int rankSelf= block->getRank(); 
	if (!(offx+offy+offz)==1) throw UbException(UB_EXARGS,"getNeighborVelocities called for diagonal directions!");
	//////get neighbor nodes, if existent
	if ((ix1==0 && offx==1) || (ix2==0 && offy==1) || (ix3==0 && offz==1))
	{
		int RankNeighborW; 
		UbTupleDouble3 orgNodeRW =  grid->getNodeCoordinates(block,  ix1, ix2, ix3);
		double xp000= ( val<1>(orgNodeRW));
		double yp000= ( val<2>(orgNodeRW));
		double zp000= ( val<3>(orgNodeRW));

		int currentLevel = block->getLevel();
		UbTupleInt3 blockIndexes = grid->getBlockIndexes(xp000,yp000, zp000,currentLevel);
		Block3DPtr blockNeighW;

		if ((val<1>(blockIndexes)!=0 && offx==1) || (val<2>(blockIndexes)!=0 && offy==1) || (val<3>(blockIndexes)!=0 && offz==1))
		{

			blockNeighW = grid->getBlock(val<1>(blockIndexes)-offx, val<2>(blockIndexes)-offy, val<3>(blockIndexes)-offz, currentLevel);

		}
		else if (offx==1 && grid->isPeriodicX1())
		{
			blockNeighW = grid->getBlock((grid->getNX1()-1), val<2>(blockIndexes), val<3>(blockIndexes), currentLevel);
		}
		else if (offy==1 && grid->isPeriodicX1())
		{
			blockNeighW = grid->getBlock(val<1>(blockIndexes),(grid->getNX2()-1),  val<3>(blockIndexes), currentLevel);
		}
		else if (offz==1 && grid->isPeriodicX1())
		{
			blockNeighW = grid->getBlock( val<1>(blockIndexes), val<2>(blockIndexes),(grid->getNX3()-1), currentLevel);
		}
		else neighNodeIsBC;

		if(blockNeighW && blockNeighW->isActive())
		{	     		     
			RankNeighborW= blockNeighW->getRank();   
		}
		else
		{

			blockNeighW = block;
			RankNeighborW= blockNeighW->getRank();
			checkInterpolation=false;
		}
		if (RankNeighborW!=rankSelf)
		{

			blockNeighW = block;
			RankNeighborW= blockNeighW->getRank();
			checkInterpolation=false;
		}

		///////////////////////////////////////
		////compute distribution at neighboring nodes from neighboring blocks

		if (checkInterpolation==false || neighNodeIsBC)
		{

			LBMKernel3DPtr kernelW = blockNeighW->getKernel();
			BCArray3D<D3Q27BoundaryCondition>& bcArrayW = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernelW->getBCProcessor())->getBCArray();          
			DistributionArray3DPtr distributionsW = kernelW->getDataSet()->getFdistributions();
			LBMReal fW2[27];
			LBMReal fW[27];
			LBMReal f0[27];
			LBMReal fE[27];
			LBMReal v0[3];
			LBMReal vW2[3];
			//distributionsW->getDistribution(fW2, std::max(ix1+2*offx,1), std::max(ix2+2*offy,1), std::max(ix3+2*offz,1));
			//distributionsW->getDistribution(fW, std::max(ix1+offx,1), std::max(ix2+offy,1), std::max(ix3+offz,1));
			//distributionsW->getDistribution(f0, std::max(ix1    ,1), std::max(ix2    ,1), std::max(ix3    ,1));
			//distributions->getDistribution(fE, std::max(ix1+offx    ,1), std::max(ix2+offy    ,1), std::max(ix3+offz    ,1)); //E:= plus 1
			distributionsW->getDistribution(fW2, std::max(ix1+2*offx,0), std::max(ix2+2*offy,0), std::max(ix3+2*offz,0));
			distributionsW->getDistribution(fW, std::max(ix1+offx,0), std::max(ix2+offy,0), std::max(ix3+offz,0));
			distributionsW->getDistribution(f0, std::max(ix1    ,0), std::max(ix2    ,0), std::max(ix3    ,0));
			distributions->getDistribution(fE, std::max(ix1+offx    ,0), std::max(ix2+offy    ,0), std::max(ix3+offz    ,0)); //E:= plus 1

			computeVelocity(fE,vE);
			computeVelocity(fW,vW);
			computeVelocity(fW2,vW2);
			computeVelocity(f0,v0);
			//second order non-symetric interpolation
			vW[0]=v0[0]*1.5-vW[0]+0.5*vW2[0];
			vW[1]=v0[1]*1.5-vW[1]+0.5*vW2[1];
			vW[2]=v0[2]*1.5-vW[2]+0.5*vW2[2];
		    //throw UbException(UB_EXARGS,"Parallel or Non-Uniform Simulation -- not yet implemented");
		}
		else
		{
			LBMKernel3DPtr kernelW = blockNeighW->getKernel();
			BCArray3D<D3Q27BoundaryCondition>& bcArrayW = boost::dynamic_pointer_cast<D3Q27ETBCProcessor>(kernelW->getBCProcessor())->getBCArray();          
			DistributionArray3DPtr distributionsW = kernelW->getDataSet()->getFdistributions();
			LBMReal fW[27];

			if (offx==1)
			{
				distributionsW->getDistribution(fW, (distributions->getNX1())-1, ix2, ix3); //moved one block backward, now get last entry
			}
			else if (offy==1)
			{
				distributionsW->getDistribution(fW, ix1,(distributions->getNX2())-1, ix3); 

			}
			else if (offz==1)
			{
				distributionsW->getDistribution(fW, ix1,ix2,distributions->getNX3()-1); 
			}
			computeVelocity(fW,vW);
		}


	}					 
	else
	{
		//data available in current block:
		LBMReal fW[27];
		distributions->getDistribution(fW, ix1-offx, ix2-offy, ix3-offz);
		computeVelocity(fW,vW);

	}
	if (checkInterpolation==true)
	{
		//in plus-direction data is available in current block because of ghost layers
		LBMReal fE[27];
		distributions->getDistribution(fE, ix1+offx, ix2+offy, ix3+offz); //E:= plus 1
		computeVelocity(fE,vE);
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void QCriterionCoProcessor::computeVelocity(LBMReal* f, LBMReal* v)
{
	//////////////////////////////////////////////////////////////////////////
	//compute x,y,z-velocity components from distribution
	//////////////////////////////////////////////////////////////////////////
	//v[xdir] = f[E] - f[W] + f[NE] - f[SW] + f[SE] - f[NW] + f[TE] - f[BW]
	//+ f[BE] - f[TW] + f[TNE] - f[TSW] + f[TSE] - f[TNW] + f[BNE] - f[BSW]
	//+ f[BSE] - f[BNW]; 

	//v[ydir] = f[N] - f[S] + f[NE] - f[SW] - f[SE] + f[NW] + f[TN] - f[BS] + f[BN]
	//- f[TS] + f[TNE] - f[TSW] - f[TSE] + f[TNW] + f[BNE] - f[BSW] - f[BSE] 
	//+ f[BNW]; 

	//v[zdir] = f[T] - f[B] + f[TE] - f[BW] - f[BE] + f[TW] + f[TN] - f[BS] - f[BN] 
	//+ f[TS] + f[TNE] + f[TSW] + f[TSE] + f[TNW] - f[BNE] - f[BSW] - f[BSE] 
	//- f[BNW];

   v[xdir] = D3Q27System::getIncompVelocityX1(f);

   v[ydir] = D3Q27System::getIncompVelocityX2(f);

   v[zdir] = D3Q27System::getIncompVelocityX3(f);
}
