#include "Simulation.h"

#include <stdio.h>
#include <vector>

#include <helper_timer.h>

#include "LBM/LB.h"
#include "Communication/Communicator.h"
#include "Communication/ExchangeData27.h"
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "GPU/GPU_Interface.h"
#include "GPU/devCheck.h"
#include "basics/utilities/UbFileOutputASCII.h"
//////////////////////////////////////////////////////////////////////////
#include "Output/MeasurePointWriter.hpp"
#include "Output/AnalysisData.hpp"
#include "Output/InterfaceDebugWriter.hpp"
#include "Output/VeloASCIIWriter.hpp"
//////////////////////////////////////////////////////////////////////////
#include "Utilities/Buffer2D.hpp"
#include "Core/StringUtilities/StringUtil.h"
//////////////////////////////////////////////////////////////////////////
#include "Init/InitLattice.h"
#include "Init/VfReader.h"
//////////////////////////////////////////////////////////////////////////
#include "FindQ/FindQ.h"
#include "FindQ/DefineBCs.h"
//////////////////////////////////////////////////////////////////////////
#include "Particles/Particles.h"
//////////////////////////////////////////////////////////////////////////
#include "Calculation/UpdateGrid27.h"
#include "Calculation/PlaneCalculations.h"
#include "Calculation/DragLift.h"
#include "Calculation/Cp.h"
#include "Calculation/Calc2ndMoments.h"
#include "Calculation/CalcMedian.h"
#include "Calculation/ForceCalculations.h"
#include "Calculation/PorousMedia.h"
//////////////////////////////////////////////////////////////////////////
#include "Restart/RestartObject.h"
//////////////////////////////////////////////////////////////////////////
#include "DataStructureInitializer/GridProvider.h"
#include "Output/DataWriter.h"
#include "Kernel/Utilities/KernelFactory/KernelFactory.h"
#include "PreProcessor/PreProcessorFactory/PreProcessorFactory.h"
#include "Kernel/Kernel.h"



std::string getFileName(const std::string& fname, int step, int myID)
{
    return std::string(fname + "_Restart_" + UbSystem::toString(myID) + "_" +  UbSystem::toString(step));
}


void Simulation::setFactories(std::shared_ptr<KernelFactory> kernelFactory, std::shared_ptr<PreProcessorFactory> preProcessorFactory)
{
	this->kernelFactory = kernelFactory;
	this->preProcessorFactory = preProcessorFactory;
}

void Simulation::addKineticEnergyAnalyzer(uint tAnalyse)
{
    this->kineticEnergyAnalyzer = std::make_shared<KineticEnergyAnalyzer>(this->para, tAnalyse);
}

void Simulation::addEnstrophyAnalyzer(uint tAnalyse)
{
    this->enstrophyAnalyzer = std::make_shared<EnstrophyAnalyzer>(this->para, tAnalyse);
}


void Simulation::init(SPtr<Parameter> para, SPtr<GridProvider> gridProvider, std::shared_ptr<DataWriter> dataWriter, std::shared_ptr<CudaMemoryManager> cudaManager)
{
   this->dataWriter = dataWriter;
   this->gridProvider = gridProvider;
   this->cudaManager = cudaManager;
   gridProvider->initalGridInformations();
   comm = vf::gpu::Communicator::getInstanz();
   this->para = para;

   devCheck(comm->mapCudaDevice(para->getMyID(), para->getNumprocs(), para->getDevices(), para->getMaxDev()));
   
   para->initLBMSimulationParameter();

   gridProvider->allocAndCopyForcing();
   gridProvider->allocAndCopyQuadricLimiters();
   if (para->getUseStreams()) {
	   gridProvider->allocArrays_fluidNodeIndices();
	   gridProvider->allocArrays_fluidNodeIndicesBorder();
   }
   gridProvider->setDimensions();
   gridProvider->setBoundingBox();

   para->setRe(para->getVelocity() * (real)1.0 / para->getViscosity());
   para->setlimitOfNodesForVTK(30000000); //max 30 Million nodes per VTK file
   if (para->getDoRestart())
       para->setStartTurn(para->getTimeDoRestart());
   else
       para->setStartTurn((unsigned int)0); //100000

   restart_object = std::make_shared<ASCIIRestartObject>();
   //////////////////////////////////////////////////////////////////////////
   output.setName(para->getFName() + StringUtil::toString<int>(para->getMyID()) + ".log");
   if(para->getMyID() == 0) output.setConsoleOut(true);
   output.clearLogFile();
   //////////////////////////////////////////////////////////////////////////
   // CUDA streams
   if (para->getUseStreams()) {
       para->getStreamManager()->launchStreams(2u);
       para->getStreamManager()->createCudaEvents();
   }
   //////////////////////////////////////////////////////////////////////////
   // 
   //output << para->getNeedInterface().at(0) << "\n";
   //output << para->getNeedInterface().at(1) << "\n";
   //output << para->getNeedInterface().at(2) << "\n";
   //output << para->getNeedInterface().at(3) << "\n";
   //output << para->getNeedInterface().at(4) << "\n";
   //output << para->getNeedInterface().at(5) << "\n";
   //////////////////////////////////////////////////////////////////////////
   //output << "      \t GridX \t GridY \t GridZ \t DistX \t DistY \t DistZ\n";
   //for (int testout=0; testout<=para->getMaxLevel();testout++)
   //{
   //   output << "Level " << testout << ":  " << para->getGridX().at(testout) << " \t " << para->getGridY().at(testout) << " \t " << para->getGridZ().at(testout) << " \t " << para->getDistX().at(testout) << " \t " << para->getDistY().at(testout) << " \t " << para->getDistZ().at(testout) << " \n";
   //}
   //////////////////////////////////////////////////////////////////////////
   output << "LB_Modell:  D3Q"<< para->getD3Qxx()          << "\n"; 
   output << "Re:         "   << para->getRe()             << "\n";
   output << "vis_ratio:  "   << para->getViscosityRatio() << "\n";
   output << "u0_ratio:   "   << para->getVelocityRatio()  << "\n";
   output << "delta_rho:  "   << para->getDensityRatio()   << "\n";
   //////////////////////////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////////////////
   cudaManager->setMemsizeGPU(0, true);
   //////////////////////////////////////////////////////////////////////////
   gridProvider->allocArrays_CoordNeighborGeo();
   gridProvider->allocArrays_BoundaryValues();
   gridProvider->allocArrays_BoundaryQs();
   gridProvider->allocArrays_OffsetScale();


   //////////////////////////////////////////////////////////////////////////
   //Kernel init
   //////////////////////////////////////////////////////////////////////////
   output << "make Kernels  " << "\n";
   kernels = kernelFactory->makeKernels(para);
   
   output << "make AD Kernels  " << "\n";
   if (para->getDiffOn())
	   adKernels = kernelFactory->makeAdvDifKernels(para);

   //////////////////////////////////////////////////////////////////////////
   //PreProcessor init
   //////////////////////////////////////////////////////////////////////////
   output << "make Preprocessors  " << "\n";
   std::vector<PreProcessorType> preProTypes = kernels.at(0)->getPreProcessorTypes();
   preProcessor = preProcessorFactory->makePreProcessor(preProTypes, para);

   //////////////////////////////////////////////////////////////////////////
   //Particles preprocessing
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalcParticle())
   {
	   rearrangeGeometry(para.get(), cudaManager.get());
	   //////////////////////////////////////////////////////////////////////////
	   allocParticles(para.get(), cudaManager.get());
	   //////////////////////////////////////////////////////////////////////////
	   ////CUDA random number generation
	   //para->cudaAllocRandomValues();

	   ////init
	   //initRandomDevice(para->getRandomState(), 
		  // para->getParD(0)->plp.numberOfParticles, 
		  // para->getParD(0)->numberofthreads);

	   ////generate random values
	   //generateRandomValuesDevice(  para->getRandomState(), 
		  // para->getParD(0)->plp.numberOfParticles, 
		  // para->getParD(0)->plp.randomLocationInit, 
		  // para->getParD(0)->numberofthreads);

	   //////////////////////////////////////////////////////////////////////////////
	   initParticles(para.get());
   }
   ////////////////////////////////////////////////////////////////////////////

   //////////////////////////////////////////////////////////////////////////
   //Allocate Memory for Drag Lift Calculation
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalcDragLift()) allocDragLift(para.get(), cudaManager.get());


   //////////////////////////////////////////////////////////////////////////
   //Allocate Memory for Plane Conc Calculation
   //////////////////////////////////////////////////////////////////////////
   //if (para->getDiffOn()) allocPlaneConc(para.get(), cudaManager.get());


   //////////////////////////////////////////////////////////////////////////
   //Median
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalcMedian())
   {
       output << "alloc Calculation for Mean Valus  " << "\n";
	   if (para->getDiffOn())	allocMedianAD(para.get(), cudaManager.get());
	   else						allocMedian(para.get(), cudaManager.get());
   }


   //////////////////////////////////////////////////////////////////////////
   //allocate memory and initialize 2nd, 3rd and higher order moments
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalc2ndOrderMoments()){  alloc2ndMoments(para.get(), cudaManager.get());         init2ndMoments(para.get());         }
   if (para->getCalc3rdOrderMoments()){  alloc3rdMoments(para.get(), cudaManager.get());         init3rdMoments(para.get());         }
   if (para->getCalcHighOrderMoments()){ allocHigherOrderMoments(para.get(), cudaManager.get()); initHigherOrderMoments(para.get()); }


   //////////////////////////////////////////////////////////////////////////
   //MeasurePoints
   //////////////////////////////////////////////////////////////////////////
   if (para->getUseMeasurePoints())
   {
	   output << "read measure points...";
	   readMeasurePoints(para.get(), cudaManager.get());
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //Porous Media
   //////////////////////////////////////////////////////////////////////////
   if (para->getSimulatePorousMedia())
   {
	   output << "define area(s) of porous media...";
	   porousMedia();
	   kernelFactory->setPorousMedia(pm);
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //enSightGold
   //////////////////////////////////////////////////////////////////////////
   //excludeGridInterfaceNodesForMirror(para, 7);
   ////output << "print case file...";
   //printCaseFile(para);
   ////output << "done.\n";
   ////output << "print geo file...";
   //printGeoFile(para, true);  //true for binary
   ////output << "done.\n";

   //////////////////////////////////////////////////////////////////////////
   //Forcing
   //////////////////////////////////////////////////////////////////////////
   ////allocVeloForForcing(para);
   //output << "new object forceCalulator  " << "\n";
   //forceCalculator = std::make_shared<ForceCalculations>(para.get());

   //////////////////////////////////////////////////////////////////////////
   //output << "define the Grid..." ;
   //defineGrid(para, comm);
   ////allocateMemory();
   //output << "done.\n";

   output << "init lattice..." ;
   initLattice(para, preProcessor, cudaManager);
   output << "done.\n";

   //output << "set geo for Q...\n" ;
   //setGeoForQ();
   //output << "done.\n";

   //if (maxlevel>1)
   //{
      //output << "find Qs...\n" ;
      //findQ27(para);
      //output << "done.\n";
   //}

   //if (para->getDiffOn()==true)
   //{
   //   output << "define TempBC...\n" ;
   //   findTempSim(para);
   //   output << "done.\n";

   //   output << "define TempVelBC...\n" ;
   //   findTempVelSim(para);
   //   output << "done.\n";

   //   output << "define TempPressBC...\n" ;
   //   findTempPressSim(para);
   //   output << "done.\n";
   //}

   //output << "find Qs-BC...\n" ;
   //findBC27(para);
   //output << "done.\n";

   //output << "find Press-BC...\n" ;
   //findPressQShip(para);
   //output << "done.\n";

   //////////////////////////////////////////////////////////////////////////
   // find indices of corner nodes for multiGPU communication
   //////////////////////////////////////////////////////////////////////////
   if (para->getDevices().size() > 2) {
       output << "Find indices of corner nodes for multiGPU communication ...";
       para->findCornerNodesCommMultiGPU();
       output << "done.\n";
   }
   //////////////////////////////////////////////////////////////////////////
   //Memory alloc for CheckPoint / Restart
   //////////////////////////////////////////////////////////////////////////
   if (para->getDoCheckPoint() || para->getDoRestart())
   {
	   output << "Alloc Memory for CheckPoint / Restart...";
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   cudaManager->cudaAllocFsForCheckPointAndRestart(lev);
	   }
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //Restart
   //////////////////////////////////////////////////////////////////////////
   if (para->getDoRestart())
   {
	   output << "Restart...\n...get the Object...\n";

		const auto name = getFileName(para->getFName(), para->getTimeDoRestart(), para->getMyID());
		restart_object->deserialize(name, para);

	   output << "...copy Memory for Restart...\n";
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   //////////////////////////////////////////////////////////////////////////
		   cudaManager->cudaCopyFsForRestart(lev);
		   //////////////////////////////////////////////////////////////////////////
		   //macroscopic values
			CalcMacSP27(para->getParD(lev)->vx_SP,       
						para->getParD(lev)->vy_SP,        
						para->getParD(lev)->vz_SP,        
						para->getParD(lev)->rho_SP, 
						para->getParD(lev)->press_SP, 
						para->getParD(lev)->geoSP,       
						para->getParD(lev)->neighborX_SP, 
						para->getParD(lev)->neighborY_SP, 
						para->getParD(lev)->neighborZ_SP,
						para->getParD(lev)->size_Mat_SP, 
						para->getParD(lev)->numberofthreads,       
						para->getParD(lev)->d0SP.f[0],    
						para->getParD(lev)->evenOrOdd);
			getLastCudaError("Kernel CalcMacSP27 execution failed"); 
			//////////////////////////////////////////////////////////////////////////
			//test...should not work...and does not
			//para->getEvenOrOdd(lev)==false;
	   }
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //Print Init
   //////////////////////////////////////////////////////////////////////////
   output << "Print files Init...";
   dataWriter->writeInit(para, cudaManager);
   if (para->getCalcParticle()) 
       copyAndPrintParticles(para.get(), cudaManager.get(), 0, true);
   output << "done.\n";

   //////////////////////////////////////////////////////////////////////////
   output << "used Device Memory: " << cudaManager->getMemsizeGPU() / 1000000.0 << " MB\n";
   //////////////////////////////////////////////////////////////////////////

   //InterfaceDebugWriter::writeInterfaceLinesDebugCF(para.get());
   //InterfaceDebugWriter::writeInterfaceLinesDebugFC(para.get());
}

void Simulation::bulk()
{

}

void Simulation::run()
{
   double ftimeE, ftimeS, fnups, durchsatz;
   float timerE, timerS;
   timerE   = 0.0f;
   timerS   = 0.0f;
   ftimeE   = 0.0f;
   ftimeS   = 0.0f;
   unsigned int t, t_prev;
   unsigned int t_MP = 0;
   //////////////////////////////////////////////////////////////////////////
   para->setStepEnsight(0);

   //turning Ship
   real Pi = (real)3.14159265358979323846;
   real delta_x_F = (real)0.1;
   real delta_t_F = (real)((double)para->getVelocity() * (double)delta_x_F / (double)3.75); 
   real delta_t_C = (real)(delta_t_F * pow(2.,para->getMaxLevel()));
   real timesteps_C = (real)(12.5 / delta_t_C);
   real AngularVelocity = (real)(12.5 / timesteps_C * Pi / 180.);
   para->setAngularVelocity(AngularVelocity);
   for (int i = 0; i<= para->getMaxLevel(); i++)
   {
	   para->getParD(i)->deltaPhi = (real)(para->getAngularVelocity()/(pow(2.,i)));
   }
   //////////////////////////////////////////////////////////////////////////

   //Timer SDK
   StopWatchInterface *sdkTimer = NULL;
   sdkCreateTimer(&sdkTimer);
   sdkStartTimer(&sdkTimer);
   //Timer Event
   cudaEvent_t start_t, stop_t;
   checkCudaErrors( cudaEventCreate(&start_t));
   checkCudaErrors( cudaEventCreate(&stop_t));
   checkCudaErrors( cudaEventRecord(start_t));

   t_prev = para->getTimeCalcMedStart();

   output << "Processing time (ms) \t Nups in Mio \t Durchsatz in GB/sec\n";

   output << "getMaxLevel = " << para->getMaxLevel() << "\n";
	////////////////////////////////////////////////////////////////////////////////
	// Time loop
	////////////////////////////////////////////////////////////////////////////////
	for(t=para->getTStart();t<=para->getTEnd();t++)
	{
        updateGrid27(para.get(), comm, cudaManager.get(), pm, 0, t, kernels);

	    ////////////////////////////////////////////////////////////////////////////////
	    //Particles
	    ////////////////////////////////////////////////////////////////////////////////
	    if (para->getCalcParticle()) propagateParticles(para.get(), t);
	    ////////////////////////////////////////////////////////////////////////////////




        ////////////////////////////////////////////////////////////////////////////////
        // run Analyzers for kinetic energy and enstrophy for TGV in 3D
        // these analyzers only work on level 0
	    ////////////////////////////////////////////////////////////////////////////////
        if (this->kineticEnergyAnalyzer || this->enstrophyAnalyzer) {
            prepareExchangeMultiGPU(para.get(), 0, -1);
            exchangeMultiGPU(para.get(), comm, cudaManager.get(), 0, -1);
        }

	    if( this->kineticEnergyAnalyzer ) this->kineticEnergyAnalyzer->run(t);
	    if( this->enstrophyAnalyzer     ) this->enstrophyAnalyzer->run(t);
	    ////////////////////////////////////////////////////////////////////////////////




        ////////////////////////////////////////////////////////////////////////////////
        //Calc Median
        ////////////////////////////////////////////////////////////////////////////////
        if (para->getCalcMedian() && ((int)t >= para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()))
        {
          for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
          {
        	  //CalcMedSP27(para->getParD(lev)->vx_SP_Med,       
        			//	  para->getParD(lev)->vy_SP_Med,        
        			//	  para->getParD(lev)->vz_SP_Med,        
        			//	  para->getParD(lev)->rho_SP_Med, 
        			//	  para->getParD(lev)->press_SP_Med, 
        			//	  para->getParD(lev)->geoSP,       
        			//	  para->getParD(lev)->neighborX_SP, 
        			//	  para->getParD(lev)->neighborY_SP, 
        			//	  para->getParD(lev)->neighborZ_SP,
        			//	  para->getParD(lev)->size_Mat_SP, 
        			//	  para->getParD(lev)->numberofthreads,       
        			//	  para->getParD(lev)->d0SP.f[0],    
        			//	  para->getParD(lev)->evenOrOdd);
        	  //getLastCudaError("CalcMacSP27 execution failed"); 
        
        	  CalcMedCompSP27(para->getParD(lev)->vx_SP_Med,       
        					  para->getParD(lev)->vy_SP_Med,        
        					  para->getParD(lev)->vz_SP_Med,        
        					  para->getParD(lev)->rho_SP_Med, 
        					  para->getParD(lev)->press_SP_Med, 
        					  para->getParD(lev)->geoSP,       
        					  para->getParD(lev)->neighborX_SP, 
        					  para->getParD(lev)->neighborY_SP, 
        					  para->getParD(lev)->neighborZ_SP,
        					  para->getParD(lev)->size_Mat_SP, 
        					  para->getParD(lev)->numberofthreads,       
        					  para->getParD(lev)->d0SP.f[0],    
        					  para->getParD(lev)->evenOrOdd);
        	  getLastCudaError("CalcMacMedCompSP27 execution failed"); 
        
          }
        }
        ////////////////////////////////////////////////////////////////////////////////




        ////////////////////////////////////////////////////////////////////////////////
        // CheckPoint
        ////////////////////////////////////////////////////////////////////////////////
        if(para->getDoCheckPoint() && para->getTimeDoCheckPoint()>0 && t%para->getTimeDoCheckPoint()==0 && t>0 && !para->overWritingRestart(t))
        {
            //////////////////////////////////////////////////////////////////////////
            //Timer SDK
            sdkStopTimer(&sdkTimer);
            sdkResetTimer(&sdkTimer);
            //////////////////////////////////////////////////////////////////////////
            //Timer Event
            checkCudaErrors( cudaEventRecord(stop_t));
            checkCudaErrors( cudaEventSynchronize(stop_t));
            
            if( para->getDoCheckPoint() )
            {
                output << "Dateien fuer CheckPoint kopieren t=" << t << "...\n";
                
                for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
                {
                    cudaManager->cudaCopyFsForCheckPoint(lev);
                }
                
                output << "Dateien fuer CheckPoint schreiben t=" << t << "...";

				const auto name = getFileName(para->getFName(), t, para->getMyID());
				restart_object->serialize(name, para);

                output << "\n fertig\n";
            }
            //////////////////////////////////////////////////////////////////////////
            //Timer SDK
            sdkStartTimer(&sdkTimer);
            //////////////////////////////////////////////////////////////////////////
            //Timer Event
            checkCudaErrors( cudaEventRecord(start_t));
        }
        //////////////////////////////////////////////////////////////////////////////





        ////////////////////////////////////////////////////////////////////////////////
        //Measure Points
        ////////////////////////////////////////////////////////////////////////////////
        //set MP-Time
        if (para->getUseMeasurePoints())
        {
            if ((t%para->getTimestepForMP()) == 0)
            {
                unsigned int valuesPerClockCycle = (unsigned int)(para->getclockCycleForMP() / para->getTimestepForMP());
                for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
                {
                    //output << "start level = " << lev << "\n";
                    LBCalcMeasurePoints27(  para->getParD(lev)->VxMP,			para->getParD(lev)->VyMP,			para->getParD(lev)->VzMP,
                    				        para->getParD(lev)->RhoMP,		    para->getParD(lev)->kMP,			para->getParD(lev)->numberOfPointskMP,
                    				        valuesPerClockCycle,				t_MP,								para->getParD(lev)->geoSP,
                    				        para->getParD(lev)->neighborX_SP,   para->getParD(lev)->neighborY_SP,	para->getParD(lev)->neighborZ_SP,
                    				        para->getParD(lev)->size_Mat_SP,	para->getParD(lev)->d0SP.f[0],		para->getParD(lev)->numberofthreads,
                    				        para->getParD(lev)->evenOrOdd);
                }
                t_MP++;
            }
            
            //Copy Measure Values
            if ((t % (unsigned int)para->getclockCycleForMP()) == 0)
            {
                for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
                {
                    cudaManager->cudaCopyMeasurePointsToHost(lev);
                    para->copyMeasurePointsArrayToVector(lev);
                    output << "\n Write MeasurePoints at level = " << lev << " and timestep = " << t << "\n";
                    for (int j = 0; j < (int)para->getParH(lev)->MP.size(); j++)
                    {
                        MeasurePointWriter::writeMeasurePoints(para.get(), lev, j, t);
                    }
                    MeasurePointWriter::calcAndWriteMeanAndFluctuations(para.get(), lev, t, para->getTStartOut());
                }
                t_MP = 0;
            }
        }
        //////////////////////////////////////////////////////////////////////////////////




        //////////////////////////////////////////////////////////////////////////////////
        ////get concentration at the plane
        //////////////////////////////////////////////////////////////////////////////////
        if (para->getDiffOn() && para->getCalcPlaneConc()) 
        {
            PlaneConcThS27( para->getParD(0)->ConcPlaneIn,
            		       para->getParD(0)->cpTopIndex,
            		       para->getParD(0)->numberOfPointsCpTop,
            		       para->getParD(0)->geoSP,       
            		       para->getParD(0)->neighborX_SP, 
            		       para->getParD(0)->neighborY_SP, 
            		       para->getParD(0)->neighborZ_SP,
            		       para->getParD(0)->size_Mat_SP, 
            		       para->getParD(0)->numberofthreads,       
            		       para->getParD(0)->d27.f[0],    
            		       para->getParD(0)->evenOrOdd);
            getLastCudaError("PlaneConcThS27 execution failed"); 
            PlaneConcThS27( para->getParD(0)->ConcPlaneOut1,
            		        para->getParD(0)->cpBottomIndex,
            		        para->getParD(0)->numberOfPointsCpBottom,
            		        para->getParD(0)->geoSP,       
            		        para->getParD(0)->neighborX_SP, 
            		        para->getParD(0)->neighborY_SP, 
            		        para->getParD(0)->neighborZ_SP,
            		        para->getParD(0)->size_Mat_SP, 
            		        para->getParD(0)->numberofthreads,       
            		        para->getParD(0)->d27.f[0],    
            		        para->getParD(0)->evenOrOdd);
            getLastCudaError("PlaneConcThS27 execution failed"); 
            PlaneConcThS27( para->getParD(0)->ConcPlaneOut2,
            		        para->getParD(0)->QPress.kN,
            		        para->getParD(0)->QPress.kQ,
            		        para->getParD(0)->geoSP,       
            		        para->getParD(0)->neighborX_SP, 
            		        para->getParD(0)->neighborY_SP, 
            		        para->getParD(0)->neighborZ_SP,
            		        para->getParD(0)->size_Mat_SP, 
            		        para->getParD(0)->numberofthreads,       
            		        para->getParD(0)->d27.f[0],    
            		        para->getParD(0)->evenOrOdd);
            getLastCudaError("PlaneConcThS27 execution failed"); 
            //////////////////////////////////////////////////////////////////////////////////
            ////Calculation of concentration at the plane
            //////////////////////////////////////////////////////////////////////////////////
            calcPlaneConc(para.get(), cudaManager.get(), 0);
        }
        //////////////////////////////////////////////////////////////////////////////////




	  ////////////////////////////////////////////////////////////////////////////////
      // File IO
      ////////////////////////////////////////////////////////////////////////////////
      //comm->startTimer();
      if(para->getTOut()>0 && t%para->getTOut()==0 && t>para->getTStartOut())
      {
		  //////////////////////////////////////////////////////////////////////////////////
		  //if (para->getParD(0)->evenOrOdd==true)  para->getParD(0)->evenOrOdd=false;
		  //else                                    para->getParD(0)->evenOrOdd=true;
		  //////////////////////////////////////////////////////////////////////////////////

		  
		 //////////////////////////////////////////////////////////////////////////
		 //Timer SDK
		 checkCudaErrors(cudaDeviceSynchronize());
		 sdkStopTimer(&sdkTimer);
		 timerS = sdkGetTimerValue(&sdkTimer);
		 sdkResetTimer(&sdkTimer);
		 ftimeS += timerS;
		 fnups = 0.0;
		 durchsatz = 0.0;
		 for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
		 {
			 fnups += 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP * pow(2.,lev) / (ftimeS*1.0E6);
			 durchsatz  +=  (27.0+1.0) * 4.0 * 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP  / (ftimeS*1.0E9);
		 }
		 output << timerS << " / " << ftimeS << " \t " <<  fnups << " \t " << durchsatz << "\n";
         //////////////////////////////////////////////////////////////////////////
		 //Timer Event
		 checkCudaErrors( cudaEventRecord(stop_t));
         checkCudaErrors( cudaEventSynchronize(stop_t));
         checkCudaErrors( cudaEventElapsedTime( &timerE, start_t, stop_t));
         ftimeE += timerE;
         fnups = 0.0;
         durchsatz = 0.0;
         for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
         {
            fnups += 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP * pow(2.,lev) / (ftimeE*1.0E6);
            durchsatz  +=  (27.0+1.0) * 4.0 * 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP  / (ftimeE*1.0E9);
         }
         output << timerE << " / " << ftimeE << " \t " <<  fnups << " \t " << durchsatz << "\n";

         if( para->getPrintFiles() )
         {
            output << "Dateien schreiben t=" << t << "...";
            for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
            {
		        //////////////////////////////////////////////////////////////////////////
		        //exchange data for valid post process
                prepareExchangeMultiGPU(para.get(), lev, -1);
		        exchangeMultiGPU(para.get(), comm, cudaManager.get(), lev, -1);
                //////////////////////////////////////////////////////////////////////////
               //if (para->getD3Qxx()==19)
               //{
                  //CalcMac(para->getParD(lev)->vx,     para->getParD(lev)->vy,       para->getParD(lev)->vz,      para->getParD(lev)->rho, 
                  //        para->getParD(lev)->geo,    para->getParD(lev)->size_Mat, para->getParD(lev)->gridNX,  para->getParD(lev)->gridNY, 
                  //        para->getParD(lev)->gridNZ, para->getParD(lev)->d0.f[0],  para->getParD(lev)->evenOrOdd);
               //}
               //else if (para->getD3Qxx()==27)
               //{
				   //if (para->getCalcMedian() && ((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()))
				   //{
					  // unsigned int tdiff = t - t_prev;
					  // CalcMacMedSP27(para->getParD(lev)->vx_SP_Med,       
				   //					  para->getParD(lev)->vy_SP_Med,        
				   //					  para->getParD(lev)->vz_SP_Med,        
				   //					  para->getParD(lev)->rho_SP_Med, 
				   //					  para->getParD(lev)->press_SP_Med, 
				   //					  para->getParD(lev)->geoSP,       
				   //					  para->getParD(lev)->neighborX_SP, 
				   //					  para->getParD(lev)->neighborY_SP, 
				   //					  para->getParD(lev)->neighborZ_SP,
				   //					  tdiff,
				   //					  para->getParD(lev)->size_Mat_SP, 
				   //					  para->getParD(lev)->numberofthreads,       
				   //					  para->getParD(lev)->evenOrOdd);
					  // getLastCudaError("CalcMacMedSP27 execution failed"); 
				   //}

				   //CalcMacSP27(para->getParD(lev)->vx_SP,       
       //                        para->getParD(lev)->vy_SP,        
       //                        para->getParD(lev)->vz_SP,        
       //                        para->getParD(lev)->rho_SP, 
       //                        para->getParD(lev)->press_SP, 
       //                        para->getParD(lev)->geoSP,       
       //                        para->getParD(lev)->neighborX_SP, 
       //                        para->getParD(lev)->neighborY_SP, 
       //                        para->getParD(lev)->neighborZ_SP,
       //                        para->getParD(lev)->size_Mat_SP, 
       //                        para->getParD(lev)->numberofthreads,       
       //                        para->getParD(lev)->d0SP.f[0],    
       //                        para->getParD(lev)->evenOrOdd);
       //            getLastCudaError("CalcMacSP27 execution failed"); 

				   
				   CalcMacCompSP27(para->getParD(lev)->vx_SP,       
								   para->getParD(lev)->vy_SP,        
								   para->getParD(lev)->vz_SP,        
								   para->getParD(lev)->rho_SP, 
								   para->getParD(lev)->press_SP, 
								   para->getParD(lev)->geoSP,       
								   para->getParD(lev)->neighborX_SP, 
								   para->getParD(lev)->neighborY_SP, 
								   para->getParD(lev)->neighborZ_SP,
								   para->getParD(lev)->size_Mat_SP, 
								   para->getParD(lev)->numberofthreads,       
								   para->getParD(lev)->d0SP.f[0],    
								   para->getParD(lev)->evenOrOdd);
                   getLastCudaError("CalcMacSP27 execution failed"); 

				   //�berschreiben mit Wandknoten
				   //SetOutputWallVelocitySP27(  para->getParD(lev)->numberofthreads,
							//				   para->getParD(lev)->vx_SP,       
							//				   para->getParD(lev)->vy_SP,        
							//				   para->getParD(lev)->vz_SP,
							//				   para->getParD(lev)->QGeom.Vx,      
							//				   para->getParD(lev)->QGeom.Vy,   
							//				   para->getParD(lev)->QGeom.Vz,
							//				   para->getParD(lev)->QGeom.kQ,      
							//				   para->getParD(lev)->QGeom.k,
							//				   para->getParD(lev)->rho_SP, 
							//				   para->getParD(lev)->press_SP, 
							//				   para->getParD(lev)->geoSP,       
							//				   para->getParD(lev)->neighborX_SP, 
							//				   para->getParD(lev)->neighborY_SP, 
							//				   para->getParD(lev)->neighborZ_SP,
							//				   para->getParD(lev)->size_Mat_SP, 
							//				   para->getParD(lev)->d0SP.f[0],    
							//				   para->getParD(lev)->evenOrOdd);
       //            getLastCudaError("SetOutputWallVelocitySP27 execution failed"); 

   				   //SetOutputWallVelocitySP27(  para->getParD(lev)->numberofthreads,
										//	   para->getParD(lev)->vx_SP,       
										//	   para->getParD(lev)->vy_SP,        
										//	   para->getParD(lev)->vz_SP,
										//	   para->getParD(lev)->Qinflow.Vx,      
										//	   para->getParD(lev)->Qinflow.Vy,   
										//	   para->getParD(lev)->Qinflow.Vz,
										//	   para->getParD(lev)->kInflowQ,      
										//	   para->getParD(lev)->Qinflow.k,
										//	   para->getParD(lev)->rho_SP, 
										//	   para->getParD(lev)->press_SP, 
										//	   para->getParD(lev)->geoSP,       
										//	   para->getParD(lev)->neighborX_SP, 
										//	   para->getParD(lev)->neighborY_SP, 
										//	   para->getParD(lev)->neighborZ_SP,
										//	   para->getParD(lev)->size_Mat_SP, 
										//	   para->getParD(lev)->d0SP.f[0],    
										//	   para->getParD(lev)->evenOrOdd);
          //         getLastCudaError("SetOutputWallVelocitySP27 execution failed"); 

				 //}

				   cudaManager->cudaCopyPrint(lev);
			   if (para->getCalcMedian())
			   {
				   cudaManager->cudaCopyMedianPrint(lev);
			   }

			   //////////////////////////////////////////////////////////////////////////
               //TODO: implement flag to write ASCII data
			   if (para->getWriteVeloASCIIfiles())
				   VeloASCIIWriter::writeVelocitiesAsTXT(para.get(), lev, t);
			   //////////////////////////////////////////////////////////////////////////
               if( this->kineticEnergyAnalyzer || this->enstrophyAnalyzer )
               {
                   std::string fname = para->getFName() + "_ID_" + StringUtil::toString<int>(para->getMyID()) + "_t_" + StringUtil::toString<int>(t);

                   if (this->kineticEnergyAnalyzer) this->kineticEnergyAnalyzer->writeToFile(fname);
                   if (this->enstrophyAnalyzer)     this->enstrophyAnalyzer->writeToFile(fname);
               }
			   //////////////////////////////////////////////////////////////////////////


			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               if (para->getDiffOn()==true)
               {
                  if (para->getDiffMod() == 7)
                  {
                     CalcMacThS7(   para->getParD(lev)->Conc, 
                                    para->getParD(lev)->geoSP,       
                                    para->getParD(lev)->neighborX_SP, 
                                    para->getParD(lev)->neighborY_SP, 
                                    para->getParD(lev)->neighborZ_SP,
                                    para->getParD(lev)->size_Mat_SP, 
                                    para->getParD(lev)->numberofthreads,       
                                    para->getParD(lev)->d7.f[0],    
                                    para->getParD(lev)->evenOrOdd);
                     getLastCudaError("CalcMacTh7 execution failed"); 
                  } 
                  else if (para->getDiffMod() == 27)
                  {
                     CalcMacThS27(  para->getParD(lev)->Conc, 
                                    para->getParD(lev)->geoSP,       
                                    para->getParD(lev)->neighborX_SP, 
                                    para->getParD(lev)->neighborY_SP, 
                                    para->getParD(lev)->neighborZ_SP,
                                    para->getParD(lev)->size_Mat_SP, 
                                    para->getParD(lev)->numberofthreads,       
                                    para->getParD(lev)->d27.f[0],    
                                    para->getParD(lev)->evenOrOdd);
                     getLastCudaError("CalcMacTh27 execution failed"); 
                  }

				  cudaManager->cudaCopyConcDH(lev);
                  //cudaMemoryCopy(para->getParH(lev)->Conc, para->getParD(lev)->Conc,  para->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost);
               }
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   ////print cp
			   //if ((para->getParH(lev)->cpTop.size() > 0) && (t > para->getTStartOut()))
			   //{
				  // printCpTopIntermediateStep(para, t, lev);
			   //}
			   ////////////////////////////////////////////////////////////////////////////////
			   //MeasurePointWriter::writeSpacialAverageForXZSlices(para, lev, t);
			   ////////////////////////////////////////////////////////////////////////////////
			   //MeasurePointWriter::writeTestAcousticXY(para, lev, t);
			   //MeasurePointWriter::writeTestAcousticYZ(para, lev, t);
			   //MeasurePointWriter::writeTestAcousticXZ(para, lev, t);
			   ////////////////////////////////////////////////////////////////////////
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////test print press mirror
			//if (t > para->getTStartOut())
			//{
			//	////////////////////////////////////////////////////////////////////////////////
			//	//Level 7
			//	CalcCPtop27(para->getParD(7)->d0SP.f[0],
			//		para->getParD(7)->cpTopIndex,
			//		para->getParD(7)->numberOfPointsCpTop,
			//		para->getParD(7)->cpPressTop,
			//		para->getParD(7)->neighborX_SP,
			//		para->getParD(7)->neighborY_SP,
			//		para->getParD(7)->neighborZ_SP,
			//		para->getParD(7)->size_Mat_SP,
			//		para->getParD(7)->evenOrOdd,
			//		para->getParD(7)->numberofthreads);
			//	//////////////////////////////////////////////////////////////////////////////////
			//	calcPressForMirror(para, 7);
			//	////////////////////////////////////////////////////////////////////////////////
			//	//Level 8
			//	CalcCPtop27(para->getParD(8)->d0SP.f[0],
			//		para->getParD(8)->cpTopIndex,
			//		para->getParD(8)->numberOfPointsCpTop,
			//		para->getParD(8)->cpPressTop,
			//		para->getParD(8)->neighborX_SP,
			//		para->getParD(8)->neighborY_SP,
			//		para->getParD(8)->neighborZ_SP,
			//		para->getParD(8)->size_Mat_SP,
			//		para->getParD(8)->evenOrOdd,
			//		para->getParD(8)->numberofthreads);
			//	//////////////////////////////////////////////////////////////////////////////////
			//	calcPressForMirror(para, 8);
			//	////////////////////////////////////////////////////////////////////////////////
			//	//print press mirror
			//	printScalars(para, false);
			//	////////////////////////////////////////////////////////////////////////////////
			//}
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//t_prev = t;

			//////////////////////////////////////////////////////////////////////////
			////Data Analysis
			////AnalysisData::writeAnalysisData(para, t);
			//AnalysisData::writeAnalysisDataX(para, t);
			//AnalysisData::writeAnalysisDataZ(para, t);
			//////////////////////////////////////////////////////////////////////////

            ////////////////////////////////////////////////////////////////////////
            //pressure difference
            ////////////////////////////////////////////////////////////////////////
			   //if (para->getMyID() == para->getPressInID())       calcPressure(para,  "in", 0);
			   //else if (para->getMyID() == para->getPressOutID()) calcPressure(para, "out", 0);
            ////////////////////////////////////////////////////////////////////////
            //flow rate
            ////////////////////////////////////////////////////////////////////////
		      //calcFlowRate(para, 0);
            ////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////
			//calculate 2nd, 3rd and higher order moments
			////////////////////////////////////////////////////////////////////////
			if (para->getCalc2ndOrderMoments())  calc2ndMoments(para.get(), cudaManager.get());
			if (para->getCalc3rdOrderMoments())  calc3rdMoments(para.get(), cudaManager.get());
			if (para->getCalcHighOrderMoments()) calcHigherOrderMoments(para.get(), cudaManager.get());
			////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////
			//calculate median on host
			////////////////////////////////////////////////////////////////////////
			if (para->getCalcMedian() && ((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()) && ((t%(unsigned int)para->getclockCycleForMP())==0))
			{
				unsigned int tdiff = t - t_prev;
				calcMedian(para.get(), tdiff);

				/////////////////////////////////
				//added for incremental averaging
				t_prev = t;
				resetMedian(para.get());
				/////////////////////////////////
			}
			////////////////////////////////////////////////////////////////////////
			dataWriter->writeTimestep(para, t);
			////////////////////////////////////////////////////////////////////////
            if (para->getCalcDragLift()) printDragLift(para.get(), cudaManager.get(), t);
			////////////////////////////////////////////////////////////////////////
			if (para->getCalcParticle()) copyAndPrintParticles(para.get(), cudaManager.get(), t, false);
			////////////////////////////////////////////////////////////////////////
			output << "done.\n";
			////////////////////////////////////////////////////////////////////////
         }
		 sdkStartTimer(&sdkTimer);
         checkCudaErrors( cudaEventRecord(start_t));
      }
	}


	//////////////////////////////////////////////////////////////////////////
	//Timer SDK
	sdkStopTimer(&sdkTimer);
	timerS = sdkGetTimerValue(&sdkTimer);
	ftimeS += timerS;
	fnups = 0.0;
	durchsatz = 0.0;
	for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	{
		fnups += 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP * pow(2.,lev) / (ftimeS*1.0E6);
		durchsatz  +=  (27.0+1.0) * 4.0 * 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP / (ftimeS*1.0E9);
	}
	output << "Processing time: " << ftimeS << "(ms)\n";
	output << "Nups in Mio: " << fnups << "\n";
	output << "Durchsatz in GB/sec: " << durchsatz << "\n";
    //////////////////////////////////////////////////////////////////////////
	//Timer Event
    checkCudaErrors( cudaEventRecord(stop_t));
    checkCudaErrors( cudaEventSynchronize(stop_t));
    checkCudaErrors( cudaEventElapsedTime( &timerE, start_t, stop_t ));
    ftimeE += timerE;
    fnups = 0.0;
    durchsatz = 0.0;
    for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
    {
       fnups += 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP * pow(2.,lev) / (ftimeE*1.0E6);
       durchsatz  +=  (27.0+1.0) * 4.0 * 1000.0 * (t-para->getTStart()) * para->getParH(lev)->size_Mat_SP / (ftimeE*1.0E9);
    }
    output << "Processing time: " << ftimeE << "(ms)\n";
    output << "Nups in Mio: " << fnups << "\n";
    output << "Durchsatz in GB/sec: " << durchsatz << "\n";
	//////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////
    // When using multiple GPUs, get Nups of all processes
	if (para->getMaxDev() > 1) {
        std::vector<double> nups = comm->gatherNUPS(fnups);
        if (comm->getPID() == 0) {
			double sum = 0;
            for (uint pid = 0; pid < nups.size(); pid++) {
                output << "Process " << pid << ": Nups in Mio: " << nups[pid] << "\n";
                sum += nups[pid];
			}
            output << "Sum of all processes: Nups in Mio: " << sum << "\n";
		}
	}
    //////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	//printDragLift(para);
	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	if (para->getDiffOn()==true) printPlaneConc(para.get(), cudaManager.get());
	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	////for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	////{
	////	if (para->getParH(lev)->cpTop.size() > 0)
	////	{
	////		printCpTop(para, lev);
	////	}
	////}
	//for (int lev = 7; lev <= 8; lev++)
	//{
	//	printCpTop(para, lev);
	//}
	////printCpTop(para);
	////printCpBottom(para);
	////printCpBottom2(para);
	////////////////////////////////////////////////////////////////////////////////

 //  //////////////////////////////////////////////////////////////////////////
 //  //Copy Measure Values
	//for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	//{
	//	output << "\n Copy MeasurePoints at level = " << lev <<"\n";
	//	para->cudaCopyMeasurePointsToHost(lev);
	//	para->copyMeasurePointsArrayToVector(lev);
	//	output << "\n Write MeasurePoints at level = " << lev <<"\n";
	//	for(int j = 0; j < (int)para->getParH(lev)->MP.size(); j++)
	//	{
	//		MeasurePointWriter::writeMeasurePoints(para, lev, j, 0);
	//	}
	//}                                                  
 //  //////////////////////////////////////////////////////////////////////////

	checkCudaErrors(cudaEventDestroy(start_t));
	checkCudaErrors(cudaEventDestroy(stop_t));
	sdkDeleteTimer(&sdkTimer);   
}

void Simulation::porousMedia()
{
	double porosity, darcySI, forchheimerSI;
	double dxLBM = 0.00390625;
	double dtLBM = 0.00000658;
	unsigned int level, geo;
	double startX, startY, startZ, endX, endY, endZ;
	//////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////
	////Test = porous media 0
	//porosity = 0.7;
	//darcySI = 137.36; //[1/s]
	//forchheimerSI = 1037.8; //[1/m]
	//level = para->getFine();
	//geo = GEO_PM_0;
	//startX = 20.0;
	//startY =  0.0;
	//startZ =  0.0;
	//endX = 40.0;
	//endY = 22.0;
	//endZ = 22.0;
	//pm[0] = new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level);
	//pm[0]->setStartCoordinates(startX, startY, startZ);
	//pm[0]->setEndCoordinates(endX, endY, endZ);
	//pm[0]->setResistanceLBM();
	//definePMarea(pm[0]);
	////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//Kondensator = porous media 0
	porosity = 0.7;
	darcySI = 137.36; //[1/s]
	forchheimerSI = 1037.8; //[1/m]
	level = para->getFine();
	geo = GEO_PM_0;
	startX = -0.715882;
	startY = -0.260942;
	startZ = -0.031321;
	endX = -0.692484;
	endY =  0.277833;
	endZ =  0.360379;
	pm.push_back(std::shared_ptr<PorousMedia>(new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level)));
	int n = (int)pm.size() - 1;
	pm.at(n)->setStartCoordinates(startX, startY, startZ);
	pm.at(n)->setEndCoordinates(endX, endY, endZ);
	pm.at(n)->setResistanceLBM();
	definePMarea(pm.at(n));
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//NT-Kuehler = porous media 1
	porosity = 0.6;
	darcySI = 149.98; //[1/s]
	forchheimerSI = 960.57; //[1/m]
	level = para->getFine();
	geo = GEO_PM_1;
	startX = -0.696146;
	startY = -0.32426;
	startZ = -0.0421345;
	endX = -0.651847;
	endY =  0.324822;
	endZ =  0.057098;
	pm.push_back(std::shared_ptr<PorousMedia>(new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level)));
	n = (int)pm.size() - 1;
	pm.at(n)->setStartCoordinates(startX, startY, startZ);
	pm.at(n)->setEndCoordinates(endX, endY, endZ);
	pm.at(n)->setResistanceLBM();
	definePMarea(pm.at(n));
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//Wasserkuehler = porous media 2
	porosity = 0.6;
	darcySI = 148.69; //[1/s]
	forchheimerSI = 629.45; //[1/m]
	level = para->getFine();
	geo = GEO_PM_2;
	startX = -0.692681;
	startY = -0.324954;
	startZ = 0.0789429;
	endX = -0.657262;
	endY =  0.32538;
	endZ =  0.400974;
	pm.push_back(std::shared_ptr<PorousMedia>(new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level)));
	n = (int)pm.size() - 1;
	pm.at(n)->setStartCoordinates(startX, startY, startZ);
	pm.at(n)->setEndCoordinates(endX, endY, endZ);
	pm.at(n)->setResistanceLBM();
	definePMarea(pm.at(n));
	//////////////////////////////////////////////////////////////////////////

}

void Simulation::definePMarea(std::shared_ptr<PorousMedia> pMedia)
{
	unsigned int counter = 0;
	unsigned int level = pMedia->getLevelPM();
	std::vector< unsigned int > nodeIDsPorousMedia;
	output << "definePMarea....find nodes \n";

	for (unsigned int i = 0; i < para->getParH(level)->size_Mat_SP; i++)
	{
		if (((para->getParH(level)->coordX_SP[i] >= pMedia->getStartX()) && (para->getParH(level)->coordX_SP[i] <= pMedia->getEndX())) &&
			((para->getParH(level)->coordY_SP[i] >= pMedia->getStartY()) && (para->getParH(level)->coordY_SP[i] <= pMedia->getEndY())) &&
			((para->getParH(level)->coordZ_SP[i] >= pMedia->getStartZ()) && (para->getParH(level)->coordZ_SP[i] <= pMedia->getEndZ())) )
		{
			if (para->getParH(level)->geoSP[i] >= GEO_FLUID)
			{
				para->getParH(level)->geoSP[i] = pMedia->getGeoID();
				nodeIDsPorousMedia.push_back(i);
				counter++;
			}
		}
	}

	output << "definePMarea....cuda copy SP \n";
	cudaManager->cudaCopySP(level);
	pMedia->setSizePM(counter);
	output << "definePMarea....cuda alloc PM \n";
	cudaManager->cudaAllocPorousMedia(pMedia.get(), level);
	unsigned int *tpmArrayIDs = pMedia->getHostNodeIDsPM();
	
	output << "definePMarea....copy vector to array \n";
	for (unsigned int j = 0; j < pMedia->getSizePM(); j++)
	{
		tpmArrayIDs[j] = nodeIDsPorousMedia[j];
	}
	
	pMedia->setHostNodeIDsPM(tpmArrayIDs);
	output << "definePMarea....cuda copy PM \n";
	cudaManager->cudaCopyPorousMedia(pMedia.get(), level);
}

void Simulation::free()
{
	// Cuda Streams
    if (para->getUseStreams()) {
        para->getStreamManager()->destroyCudaEvents();
        para->getStreamManager()->terminateStreams();
	}

	//CudaFreeHostMemory
    for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
	{
		//para->cudaFreeFull(lev);
		cudaManager->cudaFreeCoord(lev);
		cudaManager->cudaFreeSP(lev);
		if (para->getCalcMedian())
		{
			cudaManager->cudaFreeMedianSP(lev);
		}
		//para->cudaFreeVeloBC(lev);
		//para->cudaFreeWallBC(lev);
		//para->cudaFreeVeloBC(lev); 
		//para->cudaFreeInlet(lev);
		//para->cudaFreeOutlet(lev);
		//para->cudaFreeGeomBC(lev);
		//para->cudaFreePress(lev);
	}
	if (para->getMaxLevel()>1)
	{
		for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
		{
			cudaManager->cudaFreeInterfaceCF(lev);
			cudaManager->cudaFreeInterfaceFC(lev);
			cudaManager->cudaFreeInterfaceOffCF(lev);
			cudaManager->cudaFreeInterfaceOffFC(lev);
			//para->cudaFreePressX1(lev);
		}
	}
	//para->cudaFreeVeloBC(0); //level = 0
	//para->cudaFreePressBC();
	//para->cudaFreeVeloPropeller(para->getFine());
	//para->cudaFreePressX0(para->getCoarse());

	//////////////////////////////////////////////////////////////////////////
	//Temp
	if (para->getDiffOn() == true)
	{
		for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
		{
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->Conc_Full));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->Conc));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->Temp.temp));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->Temp.k));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempVel.temp));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempVel.velo));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempVel.k));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempPress.temp));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempPress.velo));
			checkCudaErrors(cudaFreeHost(para->getParH(lev)->TempPress.k));
		}
	}
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//free second order moments
	if (para->getCalc2ndOrderMoments())
	{
		for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
		{
			cudaManager->cudaFree2ndMoments(lev);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//free third order moments
	if (para->getCalc3rdOrderMoments())
	{
		for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
		{
			cudaManager->cudaFree3rdMoments(lev);
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//free higher order moments
	if (para->getCalcHighOrderMoments())
	{
		for (int lev = para->getCoarse(); lev <= para->getFine(); lev++)
		{
			cudaManager->cudaFreeHigherMoments(lev);
		}
	}
	//////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////
	//Multi GPU
	//////////////////////////////////////////////////////////////////////////
	////1D domain decomposition
	//if (para->getNumprocs() > 1)
	//{
	// for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
	// {
	//  for (unsigned int i=0; i < para->getNumberOfProcessNeighbors(lev, "send"); i++)
	//  {
	//   para->cudaFreeProcessNeighbor(lev, i);
	//  }
	// }
	//}
	//////////////////////////////////////////////////////////////////////////
	//3D domain decomposition
	if (para->getNumprocs() > 1)
	{
		for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
		{
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsX(lev, "send"); i++)
			{
				cudaManager->cudaFreeProcessNeighborX(lev, i);
			}
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsY(lev, "send"); i++)
			{
				cudaManager->cudaFreeProcessNeighborY(lev, i);
			}
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int i = 0; i < para->getNumberOfProcessNeighborsZ(lev, "send"); i++)
			{
				cudaManager->cudaFreeProcessNeighborZ(lev, i);
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	//Normals
	if (para->getIsGeoNormal()) {
		for (int lev = para->getCoarse(); lev < para->getFine(); lev++)
		{
			cudaManager->cudaFreeGeomNormals(lev);
		}
	}
	//////////////////////////////////////////////////////////////////////////

    delete comm;

}