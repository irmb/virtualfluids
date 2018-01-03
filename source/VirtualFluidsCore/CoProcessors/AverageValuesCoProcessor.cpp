#include "AverageValuesCoProcessor.h"

#include "LBMKernel.h"
#include "BCProcessor.h"

#include "basics/writer/WbWriterVtkXmlASCII.h"

#include "DataSet3D.h"
#include "WbWriter.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "UbScheduler.h"
#include "Communicator.h"
#include "BCArray3D.h"

using namespace std;

AverageValuesCoProcessor::AverageValuesCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
AverageValuesCoProcessor::AverageValuesCoProcessor(Grid3DPtr grid, const std::string& path,	WbWriter* const writer, 
   UbSchedulerPtr s, UbSchedulerPtr Avs, UbSchedulerPtr rsMeans, UbSchedulerPtr rsRMS, bool restart)
	                                                   : CoProcessor(grid, s),
	                                                   averageScheduler(Avs),
	                                                   resetSchedulerMeans(rsMeans),
	                                                   resetSchedulerRMS(rsRMS),
	                                                   path(path),
	                                                   writer(writer)
{
   resetStepMeans = (int)rsMeans->getMinBegin();
   resetStepRMS = (int)rsRMS->getMinBegin();
   averageInterval = (double)Avs->getMinStep();

   gridRank  = grid->getRank();
   minInitLevel = this->grid->getCoarsestInitializedLevel();
   maxInitLevel = this->grid->getFinestInitializedLevel();

   blockVector.resize(maxInitLevel+1);

   for(int level = minInitLevel; level<=maxInitLevel;level++)
   {
      grid->getBlocks(level, gridRank, true, blockVector[level]);
      
      if (blockVector[level].size() > 0)
         compressible = blockVector[level][0]->getKernel()->getCompressible();

      if (!restart)
      {
         for(Block3DPtr block : blockVector[level])
         {
            UbTupleInt3 nx = grid->getBlockNX();
            AverageValuesArray3DPtr averageValues = AverageValuesArray3DPtr(new AverageValuesArray3D(11, val<1>(nx)+1, val<2>(nx)+1, val<3>(nx)+1, 0.0));
            block->getKernel()->getDataSet()->setAverageValues(averageValues);
         }
      }
   }

   // for musis special use
	//initPlotDataZ(0.0);
	//restartStep = 0.0;
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::process(double step)
{
	//resetRMS(step);
	if(resetSchedulerRMS->isDue(step) )
		resetDataRMS(step);

	//reset(step);
	if(resetSchedulerMeans->isDue(step) )
		resetDataMeans(step);

	if(averageScheduler->isDue(step) ){
		calculateAverageValues(step);
			// for musis special use
			//collectPlotDataZ(step);
	}
	if(scheduler->isDue(step) ){
			collectData(step);

		}

		UBLOG(logDEBUG3, "AverageValuesCoProcessor::update:" << step);
}

void AverageValuesCoProcessor::resetDataRMS(double step)
{
	resetStepRMS=(int)step;

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		for(Block3DPtr block : blockVector[level])
		{
			if (block)
			{
				ILBMKernelPtr kernel = block->getKernel();
				BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();          
				DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
				AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();

				int minX1 = 0;
				int minX2 = 0;
				int minX3 = 0;

				int maxX1 = int(distributions->getNX1());
				int maxX2 = int(distributions->getNX2());
				int maxX3 = int(distributions->getNX3());

				for(int ix3=minX3; ix3<maxX3-1; ix3++)
				{
					for(int ix2=minX2; ix2<maxX2-1; ix2++)
					{
						for(int ix1=minX1; ix1<maxX1-1; ix1++)
						{
							if(!bcArray->isUndefined(ix1,ix2,ix3) && !bcArray->isSolid(ix1,ix2,ix3))
							{
								//////////////////////////////////////////////////////////////////////////
								//compute average values
								//////////////////////////////////////////////////////////////////////////
								(*av)(AvVxx,ix1,ix2,ix3) = 0.0;
								(*av)(AvVyy,ix1,ix2,ix3) = 0.0;
								(*av)(AvVzz,ix1,ix2,ix3) = 0.0;
                        (*av)(AvVxy,ix1,ix2,ix3) = 0.0;
                        (*av)(AvVxz,ix1,ix2,ix3) = 0.0;
                        (*av)(AvVyz,ix1,ix2,ix3) = 0.0;
                        (*av)(AvPrms,ix1,ix2,ix3) = 0.0;
								//////////////////////////////////////////////////////////////////////////
							}
						}
					}
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::resetDataMeans(double step)
{
	resetStepMeans=(int)step;

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		for(Block3DPtr block : blockVector[level])
		{
			if (block)
			{
				ILBMKernelPtr kernel = block->getKernel();
				BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();          
				DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
				AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();

				int minX1 = 0;
				int minX2 = 0;
				int minX3 = 0;

				int maxX1 = int(distributions->getNX1());
				int maxX2 = int(distributions->getNX2());
				int maxX3 = int(distributions->getNX3());

				for(int ix3=minX3; ix3<maxX3-1; ix3++)
				{
					for(int ix2=minX2; ix2<maxX2-1; ix2++)
					{
						for(int ix1=minX1; ix1<maxX1-1; ix1++)
						{
							if(!bcArray->isUndefined(ix1,ix2,ix3) && !bcArray->isSolid(ix1,ix2,ix3))
							{
								//////////////////////////////////////////////////////////////////////////
								//compute average values
								//////////////////////////////////////////////////////////////////////////
								(*av)(AvVx,ix1,ix2,ix3) = 0.0;
								(*av)(AvVy,ix1,ix2,ix3) = 0.0;
								(*av)(AvVz,ix1,ix2,ix3) = 0.0;
                        (*av)(AvP,ix1,ix2,ix3) = 0.0;
								//////////////////////////////////////////////////////////////////////////
							}
						}
					}
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::collectData(double step)
{
	int istep = int(step);

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		for(Block3DPtr block : blockVector[level])
		{
			if (block)
			{
				addData(block);
			}
		}
	}

   string pfilePath, partPath, subfolder, cfilePath;
   subfolder = "av"+UbSystem::toString(istep);
   pfilePath =  path+"/av/"+subfolder;
   cfilePath =  path+"/av/av_collection";
   partPath = pfilePath+"/av"+UbSystem::toString(gridRank)+ "_" + UbSystem::toString(istep);

   string partName = writer->writeOctsWithNodeData(partPath,nodes,cells,datanames,data);
   size_t found=partName.find_last_of("/");
   string piece = partName.substr(found+1);
   piece = subfolder + "/" + piece;

   vector<string> cellDataNames;
   CommunicatorPtr comm = Communicator::getInstance();
   vector<string> pieces = comm->gather(piece);
   if (comm->getProcessID() == comm->getRoot())
   {
      string pname = WbWriterVtkXmlASCII::getInstance()->writeParallelFile(pfilePath,pieces,datanames,cellDataNames);
      found=pname.find_last_of("/");
      piece = pname.substr(found+1);

      vector<string> filenames;
      filenames.push_back(piece);
      if (step == CoProcessor::scheduler->getMinBegin())
      {
         WbWriterVtkXmlASCII::getInstance()->writeCollection(cfilePath,filenames,istep,false);
      } 
      else
      {
         WbWriterVtkXmlASCII::getInstance()->addFilesToCollection(cfilePath,filenames,istep,false);
      }
      UBLOG(logINFO,"AverageValuesCoProcessor step: " << istep);
   }

	clearData();
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::clearData()
{
	nodes.clear();
	cells.clear();
	datanames.clear();
	data.clear();
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::addData(const Block3DPtr block)
{
	UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
	UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
	UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
	double         dx           = grid->getDeltaX(block);

	//Diese Daten werden geschrieben:
	datanames.resize(0);
	datanames.push_back("AvVx");
   datanames.push_back("AvVy");
   datanames.push_back("AvVz");
	datanames.push_back("AvVxx");
	datanames.push_back("AvVyy");
	datanames.push_back("AvVzz");
   datanames.push_back("AvVxy");
   datanames.push_back("AvVxz");
   datanames.push_back("AvVyz");
   datanames.push_back("AvP");
   datanames.push_back("AvPrms");


	data.resize(datanames.size());

	ILBMKernelPtr kernel = block->getKernel();
	BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();          
	DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
	AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();
	//int ghostLayerWidth = kernel->getGhostLayerWidth();

	//knotennummerierung faengt immer bei 0 an!
	int SWB,SEB,NEB,NWB,SWT,SET,NET,NWT;

	int minX1 = 0;
	int minX2 = 0;
	int minX3 = 0;

	int maxX1 = int(distributions->getNX1());
	int maxX2 = int(distributions->getNX2());
	int maxX3 = int(distributions->getNX3());

	//nummern vergeben und node vector erstellen + daten sammeln
	CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3,-1);

	maxX1 -= 2;
	maxX2 -= 2;
	maxX3 -= 2;

	//D3Q27BoundaryConditionPtr bcPtr;

	int nr = (int)nodes.size();

	for(int ix3=minX3; ix3<=maxX3; ix3++)
	{
		for(int ix2=minX2; ix2<=maxX2; ix2++)
		{
			for(int ix1=minX1; ix1<=maxX1; ix1++)
			{
				if(!bcArray->isUndefined(ix1,ix2,ix3) && !bcArray->isSolid(ix1,ix2,ix3))
				{
					int index = 0;
					nodeNumbers(ix1,ix2,ix3) = nr++;
					nodes.push_back( makeUbTuple(float(val<1>(org) - val<1>(nodeOffset) + ix1*dx),
						float(val<2>(org) - val<2>(nodeOffset) + ix2*dx),
						float(val<3>(org) - val<3>(nodeOffset) + ix3*dx)) );

					LBMReal vx=(*av)(AvVx,ix1,ix2,ix3);
					LBMReal vy=(*av)(AvVy,ix1,ix2,ix3);
					LBMReal vz=(*av)(AvVz,ix1,ix2,ix3);
               
               LBMReal vxx=(*av)(AvVxx,ix1,ix2,ix3);
               LBMReal vyy=(*av)(AvVyy,ix1,ix2,ix3);
               LBMReal vzz=(*av)(AvVzz,ix1,ix2,ix3);
               
               LBMReal vxy=(*av)(AvVxy,ix1,ix2,ix3);
               LBMReal vxz=(*av)(AvVxz,ix1,ix2,ix3);
               LBMReal vyz=(*av)(AvVyz,ix1,ix2,ix3);

               LBMReal vp=(*av)(AvP,ix1,ix2,ix3);
               LBMReal vprms=(*av)(AvPrms,ix1,ix2,ix3);
 
					
					data[index++].push_back(vx);
               data[index++].push_back(vy);
               data[index++].push_back(vz);

					data[index++].push_back(vxx);
					data[index++].push_back(vyy);
					data[index++].push_back(vzz);

               data[index++].push_back(vxy);
               data[index++].push_back(vxz);
               data[index++].push_back(vyz);

               data[index++].push_back(vp);
               data[index++].push_back(vprms);
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
					cells.push_back( makeUbTuple(SWB,SEB,NEB,NWB,SWT,SET,NET,NWT) );
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////
void AverageValuesCoProcessor::calculateAverageValues(double timeStep)
{
	using namespace D3Q27System;

   //Funktionszeiger
   calcMacros = NULL;
   if (compressible)
   {
      calcMacros = &calcCompMacroscopicValues;
   }
   else
   {
      calcMacros = &calcIncompMacroscopicValues;
   }

	LBMReal f[27];

	for(int level = minInitLevel; level<=maxInitLevel;level++)
	{
		for(Block3DPtr block : blockVector[level])
		{
			if (block)
			{
				ILBMKernelPtr kernel = block->getKernel();
				BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();          
				DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions(); 
				AverageValuesArray3DPtr av = kernel->getDataSet()->getAverageValues();

				int minX1 = 0;
				int minX2 = 0;
				int minX3 = 0;

				int maxX1 = int(distributions->getNX1());
				int maxX2 = int(distributions->getNX2());
				int maxX3 = int(distributions->getNX3());

				maxX1 -= 2;
				maxX2 -= 2;
				maxX3 -= 2;

				for(int ix3=minX3; ix3<=maxX3; ix3++)
				{
					for(int ix2=minX2; ix2<=maxX2; ix2++)
					{
						for(int ix1=minX1; ix1<=maxX1; ix1++)
						{
							if(!bcArray->isUndefined(ix1,ix2,ix3) && !bcArray->isSolid(ix1,ix2,ix3))
							{
								//////////////////////////////////////////////////////////////////////////
								//read distribution
								////////////////////////////////////////////////////////////////////////////
								distributions->getDistribution(f, ix1, ix2, ix3);
								//////////////////////////////////////////////////////////////////////////
								//compute velocity
								//////////////////////////////////////////////////////////////////////////
                        LBMReal vx,vy,vz,rho;
                        calcMacros(f,rho,vx,vy,vz);
                        double press = D3Q27System::calcPress(f,rho,vx,vy,vz);

								//////////////////////////////////////////////////////////////////////////
								//compute average values
								//////////////////////////////////////////////////////////////////////////

								LBMReal timeStepAfterResetRMS=(double)(timeStep-resetStepRMS)/((double)averageInterval);
								LBMReal timeStepAfterResetMeans=(double)(timeStep-resetStepMeans)/((double)averageInterval);

                        //mean velocity
                        (*av)(AvVx,ix1,ix2,ix3) = ((*av)(AvVx,ix1,ix2,ix3)*timeStepAfterResetMeans + vx)/(timeStepAfterResetMeans+1.0);
                        (*av)(AvVy,ix1,ix2,ix3) = ((*av)(AvVy,ix1,ix2,ix3)*timeStepAfterResetMeans + vy)/(timeStepAfterResetMeans+1.0);
                        (*av)(AvVz,ix1,ix2,ix3) = ((*av)(AvVz,ix1,ix2,ix3)*timeStepAfterResetMeans + vz)/(timeStepAfterResetMeans+1.0);

                        //rms
								(*av)(AvVxx,ix1,ix2,ix3) = ((vx-(*av)(AvVx,ix1,ix2,ix3))*(vx-(*av)(AvVx,ix1,ix2,ix3)) +
									(*av)(AvVxx,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);
								(*av)(AvVyy,ix1,ix2,ix3) = ((vy-(*av)(AvVy,ix1,ix2,ix3))*(vy-(*av)(AvVy,ix1,ix2,ix3)) +
									(*av)(AvVyy,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);
								(*av)(AvVzz,ix1,ix2,ix3) = ((vz-(*av)(AvVz,ix1,ix2,ix3))*(vz-(*av)(AvVz,ix1,ix2,ix3)) +
									(*av)(AvVzz,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);

                        //cross-correlations
                        (*av)(AvVxy,ix1,ix2,ix3) = ((vx-(*av)(AvVx,ix1,ix2,ix3))*(vy-(*av)(AvVy,ix1,ix2,ix3)) +
                           (*av)(AvVxy,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);
                        (*av)(AvVxz,ix1,ix2,ix3) = ((vx-(*av)(AvVx,ix1,ix2,ix3))*(vz-(*av)(AvVz,ix1,ix2,ix3)) +
                           (*av)(AvVxz,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);
                        (*av)(AvVyz,ix1,ix2,ix3) = ((vy-(*av)(AvVy,ix1,ix2,ix3))*(vz-(*av)(AvVz,ix1,ix2,ix3)) +
                           (*av)(AvVyz,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);

                        //mean and rms press
                        (*av)(AvP,ix1,ix2,ix3) = ((*av)(AvP,ix1,ix2,ix3)*timeStepAfterResetMeans + press)/(timeStepAfterResetMeans+1.0);
                        (*av)(AvPrms,ix1,ix2,ix3) = ((press-(*av)(AvP,ix1,ix2,ix3))*(press-(*av)(AvP,ix1,ix2,ix3)) +
                           (*av)(AvPrms,ix1,ix2,ix3)*timeStepAfterResetRMS)/(timeStepAfterResetRMS+1.0);

								//////////////////////////////////////////////////////////////////////////
							}
						}
					}
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////
//void AverageValuesCoProcessor::initPlotData(double step)
//{
//   CommunicatorPtr comm = Communicator::getInstance();
//	if (comm->getProcessID() == comm->getRoot())
//	{
//		std::ofstream ostr;
//		string fname = path + "_PlotData_" + UbSystem::toString(step) + ".txt"; 
//		ostr.open(fname.c_str(), std::ios_base::out);
//		if(!ostr)
//		{ 
//			ostr.clear();
//			string path = UbSystem::getPathFromString(fname);
//			if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out);}
//			if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
//		}
//		ostr << "Time"<< "\t" <<"Ref.Time"<<"\t"<< "Z_Coor"<< "\t" << "Pore fraction" << "\t";
//		ostr << "Vx"  << "\t" << "Vy" << "\t" << "Vz" << "\t";
//		ostr << "TSx" << "\t" << "TSy"<< "\t" << "TSz"<< "TSxz";
//		ostr << endl;
//		ostr.close();
//	}
//}
//////////////////////////////////////////////////////////////////////////////
//void AverageValuesCoProcessor::collectPlotData(double step)
//{
//
//	double hminX1 = 0.9;
//	double hminX2 = 0.0;
//	double hmaxX1 = 0.95;
//	double hmaxX2 = 0.01; //systemabmessungen world units
//
//	// 3 level platte standard:
//	double hX3_level[] = {0.305, 0.309,0.3365,0.35};
//	//0.004, 0,0365,0.045
//	//musis: 3 level refinement
//	//double hX3_level[] = {0.42, 0.28, 0.105, 0.0}; //refinement coords
//	                    //bsislevel von 0.42-0.28,... (level 0 bis 2 , 3 insgesamt)
//	//musis: 4 level refinement
//	//double hX3_level[] = {0.42, 0.3, 0.195, 0.078, 0.0};
//	//musis: 5 level refinement
//	//double hX3_level[] = {0.396, 0.28, 0.18, 0.08, 0.006, 0.0};
//
//	ostringstream Str;
//	Str << step;
//	string step2string(Str.str());
//	string fname = path + "_PlotZ_" + step2string + ".txt"; 
//
//
//	for(int level = minInitLevel; level<=maxInitLevel;level++)
//	{
//		double dx = grid->getDeltaX(level);
//
//		for (double hi =hX3_level[level]; hi >= hX3_level[level+1]; hi=hi-dx ){
//			D3Q27IntegrateValuesHelper h1(grid, comm, 
//				hminX1, hminX2, hi, 
//				hmaxX1, hmaxX2, hi-dx);
//
//			h1.calculateAV();
//			double nn1 = h1.getNumberOfNodes();
//			double ns1 = h1.getNumberOfSolids();
//			if (nn1 > 0.0){
//				// get data and write into txt files
//				if (comm->getProcessID() == comm->getRoot())
//				{
//					int istep = static_cast<int>(step);
//					std::ofstream ostr;
//
//					double AvVx1 = h1.getAvVx1()/nn1;
//					double AvVx2 = h1.getAvVx2()/nn1;
//					double AvVx3 = h1.getAvVx3()/nn1;
//
//					double AvTSx1 = h1.getTSx1()/nn1;
//					double AvTSx2 = h1.getTSx2()/nn1;
//					double AvTSx3 = h1.getTSx3()/nn1;
//
//					double AvTSx1x3 = h1.getTSx1x3()/nn1;
//
//					ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
//					if(!ostr)
//					{ 
//						ostr.clear();
//						string path = UbSystem::getPathFromString(fname);
//						if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
//						if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
//					}
//					ostr << istep << "\t" << resetStep << "\t" << hi+0.5*dx << "\t" << nn1/(nn1+ns1)*100.0 << "%\t";
//					ostr << AvVx1 << "\t" << AvVx2 << "\t" << AvVx3 << "\t";
//					ostr << AvTSx1<< "\t" << AvTSx2<< "\t" << AvTSx3<< "\t" << AvTSx1x3;
//					ostr << endl;
//					ostr.close();
//
//				}
//			}
//		}
//
//	}
//}

