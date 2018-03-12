#include "Simulation.h"
#include "LBM/LB.h"
#include "Communication/Communicator.h"
#include "Communication/ExchangeData27.h"
#include "Parameter/Parameter.h"
#include "GPU/GPU_Interface.h"
#include "GPU/devCheck.h"
#include "basics/utilities/UbFileOutputASCII.h"
//////////////////////////////////////////////////////////////////////////
#include "Input/ConfigFile.h"
#include "Input/VtkXmlReader.hpp"
#include "Input/VtkGeometryReader.h"
#include "Input/kFullReader.h"
#include "Input/PositionReader.h"
//////////////////////////////////////////////////////////////////////////
#include "Output/WriteData.h"
#include "Output/MeasurePointWriter.hpp"
#include "Output/AnalysisData.hpp"
//////////////////////////////////////////////////////////////////////////
#include "Utilities/Buffer2D.hpp"
#include "Utilities/StringUtil.hpp"
//////////////////////////////////////////////////////////////////////////
#include "Init/InitLattice.h"
#include "Init/DefineGrid.h"
#include "Init/SetParameter.h"
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
//CUDA
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>
////random numbers
//#include <curand.h>
//#include <curand_kernel.h>
//////////////////////////////////////////////////////////////////////////
#include <boost/foreach.hpp>
#include <stdio.h>
#include <vector>
#include "Restart/RestartPostprocessor.h"
//////////////////////////////////////////////////////////////////////////
#include "DataStructureInitializer/GridProvider.h"

Simulation::Simulation()
{

}

Simulation::~Simulation()
{

}

void Simulation::init(SPtr<Parameter> para, SPtr<GridProvider> gridProvider)
{
   this->gridProvider = gridProvider;
   gridProvider->initalGridInformations();
   comm = SPtr<Communicator>(Communicator::getInstanz());
   this->para = para;

   para->setMyID(comm->getPID());
   para->setNumprocs(comm->getNummberOfProcess());
   devCheck(comm->mapCudaDevice(para->getMyID(), para->getNumprocs(), para->getDevices(), para->getMaxDev()));

   gridProvider->allocAndCopyForcing();
   gridProvider->setDimensions();
   gridProvider->setBoundingBox();

   para->initParameter();
   para->setRe(para->getVelocity() * (real)1.0 / para->getViscosity());
   para->setPhi((real) 0.0);
   para->setlimitOfNodesForVTK(30000000); //max 30 Million nodes per VTK file
   if (para->getDoRestart())
       para->setStartTurn(para->getTimeDoRestart());
   else
       para->setStartTurn((unsigned int)0); //100000

   //Restart object
   restObj = new RestartObject();
   rest = new RestartPostprocessor(restObj, RestartPostprocessor::TXT);
   //////////////////////////////////////////////////////////////////////////
   output.setName(para->getFName() + StringUtil::toString<int>(para->getMyID()) + ".log");
   if(para->getMyID() == 0) output.setConsoleOut(true);
   output.clearLogFile();
   //////////////////////////////////////////////////////////////////////////
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
   para->setMemsizeGPU(0, true);
   //////////////////////////////////////////////////////////////////////////
   gridProvider->allocArrays_CoordNeighborGeo();
   gridProvider->allocArrays_BoundaryValues();
   gridProvider->allocArrays_BoundaryQs();
   gridProvider->allocArrays_OffsetScale();


   //////////////////////////////////////////////////////////////////////////
   //Particles preprocessing
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalcParticle())
   {
	   rearrangeGeometry(para.get());
	   //////////////////////////////////////////////////////////////////////////
	   allocParticles(para.get());
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
   //allocDragLift(para);


   //////////////////////////////////////////////////////////////////////////
   //Allocate Memory for Plane Conc Calculation
   //////////////////////////////////////////////////////////////////////////
   if (para->getDiffOn()) allocPlaneConc(para.get());


   //////////////////////////////////////////////////////////////////////////
   //Median
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalcMedian())   allocMedian(para.get());


   //////////////////////////////////////////////////////////////////////////
   //allocate memory and initialize 2nd, 3rd and higher order moments
   //////////////////////////////////////////////////////////////////////////
   if (para->getCalc2ndOrderMoments()){  alloc2ndMoments(para.get());         init2ndMoments(para.get());         }
   if (para->getCalc3rdOrderMoments()){  alloc3rdMoments(para.get());         init3rdMoments(para.get());         }
   if (para->getCalcHighOrderMoments()){ allocHigherOrderMoments(para.get()); initHigherOrderMoments(para.get()); }


   //////////////////////////////////////////////////////////////////////////
   //MeasurePoints
   //////////////////////////////////////////////////////////////////////////
   if (para->getUseMeasurePoints())
   {
	   output << "read measure points...";
	   readMeasurePoints(para.get());
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //Porous Media
   //////////////////////////////////////////////////////////////////////////
   if (para->getSimulatePorousMedia())
   {
	   output << "define area(s) of porous media...";
	   porousMedia();
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
   //allocVeloForForcing(para);
   forceCalculator = new ForceCalculations(para.get());

   //////////////////////////////////////////////////////////////////////////
   //output << "define the Grid..." ;
   //defineGrid(para, comm);
   ////allocateMemory();
   //output << "done.\n";

   output << "init lattice..." ;
   initLattice(para);
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
   //Memory alloc for CheckPoint / Restart
   //////////////////////////////////////////////////////////////////////////
   if (para->getDoCheckPoint() || para->getDoRestart())
   {
	   output << "Alloc Memory for CheckPoint / Restart...";
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   para->cudaAllocFsForCheckPointAndRestart(lev);
	   }
	   output << "done.\n";
   }

   //////////////////////////////////////////////////////////////////////////
   //Restart
   //////////////////////////////////////////////////////////////////////////
   if (para->getDoRestart())
   {
	   output << "Restart...\n...get the Object...\n";
	   restObj = rest->restart(para->getFName(), para->getTimeDoRestart(), para->getMyID());
	   output << "...load...\n";
	   restObj->load(para.get());
	   //para = rest->restart(para->getTimeDoRestart());
	   output << "...copy Memory for Restart...\n";
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   //////////////////////////////////////////////////////////////////////////
		   para->cudaCopyFsForRestart(lev);
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
   writeInit(para);
   if (para->getCalcParticle()) 
       copyAndPrintParticles(para.get(), 0, true);
   output << "done.\n";

   //////////////////////////////////////////////////////////////////////////
   output << "used Device Memory: " << para->getMemsizeGPU() / 1000000.0 << " MB\n";
   //////////////////////////////////////////////////////////////////////////
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
   //////////////////////////////////////////////////////////////////////////
   real visSponge = 0.001f;
   real omegaSponge = 1.f / ( 3.f * visSponge + 0.5f );
   //////////////////////////////////////////////////////////////////////////
   //turning Ship
   real Pi = (real)3.14159265358979323846;
   real delta_x_F = (real)0.1;
   real delta_t_F = (real)(para->getVelocity() * delta_x_F / 3.75); 
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
		getLastCudaError("before starting a kernel we get an execution failed");
		if (para->getMaxLevel()>=1)
         {
            updateGrid27(para.get(), comm.get(), pm, 1, para->getMaxLevel(), t);
         }
         ////////////////////////////////////////////////////////////////////////////////
         // Collision and Propagation
         ////////////////////////////////////////////////////////////////////////////////      
		 //if (t>para->getStartTurn()){
			// //////////////////////////////////////////////////////////////////////////
			// KernelKum1hSP27(    para->getParD(0)->numberofthreads,       
			//					 para->getParD(0)->omega,
			//					 para->getParD(0)->deltaPhi,
			//					 para->getAngularVelocity(),
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->coordX_SP, 
			//					 para->getParD(0)->coordY_SP, 
			//					 para->getParD(0)->coordZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->size_Mat_SP,  
			//					 para->getParD(0)->evenOrOdd); 
			// getLastCudaError("KernelCasSPKum27 execution failed");
			// //////////////////////////////////////////////////////////////////////////
			// QVelDevice1h27(    para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
			//					para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
			//					para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
			//					para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,          
			//					para->getPhi(),                    para->getAngularVelocity(),
			//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//					para->getParD(0)->coordX_SP,       para->getParD(0)->coordY_SP,    para->getParD(0)->coordZ_SP,
			//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		 //     getLastCudaError("QVelDev27 execution failed");
			//  //////////////////////////////////////////////////////////////////////////
		 //}
		 //else{
			// //////////////////////////////////////////////////////////////////////////
			// KernelKum1hSP27(    para->getParD(0)->numberofthreads,       
			//					 para->getParD(0)->omega,
			//					 (real)0.0,
			//					 (real)0.0,
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->coordX_SP, 
			//					 para->getParD(0)->coordY_SP, 
			//					 para->getParD(0)->coordZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->size_Mat_SP,  
			//					 para->getParD(0)->evenOrOdd); 
			// getLastCudaError("KernelCasSPKum27 execution failed");
			// //////////////////////////////////////////////////////////////////////////
			// QVelDevice1h27(    para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
			//					para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
			//					para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
			//					para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,          
			//					para->getPhi(),                    (real)0.0,
			//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//					para->getParD(0)->coordX_SP,       para->getParD(0)->coordY_SP,    para->getParD(0)->coordZ_SP,
			//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		 //     getLastCudaError("QVelDev27 execution failed");
			//  //////////////////////////////////////////////////////////////////////////

		 //}
		 //////////////////////////////////////////////////////////////////////////
		 //All 4
		 //KernelCumulantD3Q27All4(para->getParD(0)->numberofthreads,
		 //						 para->getParD(0)->omega, 
		 //						 para->getParD(0)->geoSP, 
		 //						 para->getParD(0)->neighborX_SP, 
		 //						 para->getParD(0)->neighborY_SP, 
		 //						 para->getParD(0)->neighborZ_SP,
		 //						 para->getParD(0)->d0SP.f[0],    
		 //						 para->getParD(0)->size_Mat_SP,
		 //						 0,
		 //						 para->getForcesDev(),
		 //						 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCumulantD3Q27All4 execution failed");
		 ////////////////////////////////////////////////////////////////////////// 
		 //F3
		 //KernelCumulantD3Q27F3( para->getParD(0)->numberofthreads,
			//					 para->getParD(0)->omega, 
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->g6.g[0],    
			//					 para->getParD(0)->size_Mat_SP,
			//					 0,
			//					 para->getForcesDev(),
			//					 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelKumAA2016CompSP27 execution failed");
		 //////////////////////////////////////////////////////////////////////////
		 //comp
		 //KernelKumAA2016CompBulkSP27(para->getParD(0)->numberofthreads,       
			//						 para->getParD(0)->omega, 
			//						 para->getParD(0)->geoSP, 
			//						 para->getParD(0)->neighborX_SP, 
			//						 para->getParD(0)->neighborY_SP, 
			//						 para->getParD(0)->neighborZ_SP,
			//						 para->getParD(0)->d0SP.f[0],    
			//						 para->getParD(0)->size_Mat_SP,
			//						 para->getParD(0)->size_Array_SP,
			//						 0,
			//						 para->getForcesDev(),
			//						 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelKumAA2016CompSP27 execution failed");
		 //KernelKumAA2016CompSP27(para->getParD(0)->numberofthreads,       
			//					 para->getParD(0)->omega, 
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->size_Mat_SP,
			//					 0,
			//					 para->getForcesDev(),
			//					 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelKumAA2016CompSP27 execution failed");
		 //KernelBGKPlusCompSP27(para->getParD(0)->numberofthreads,       
		 //					   para->getParD(0)->omega, 
		 //					   para->getParD(0)->geoSP, 
		 //					   para->getParD(0)->neighborX_SP, 
		 //					   para->getParD(0)->neighborY_SP, 
		 //					   para->getParD(0)->neighborZ_SP,
		 //					   para->getParD(0)->d0SP.f[0],    
		 //					   para->getParD(0)->size_Mat_SP,  
		 //					   para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelBGKPlusSP27 execution failed");
		 //KernelBGKCompSP27(para->getParD(0)->numberofthreads,       
		 //				   para->getParD(0)->omega, 
		 //				   para->getParD(0)->geoSP, 
		 //				   para->getParD(0)->neighborX_SP, 
		 //				   para->getParD(0)->neighborY_SP, 
		 //				   para->getParD(0)->neighborZ_SP,
		 //				   para->getParD(0)->d0SP.f[0],    
		 //				   para->getParD(0)->size_Mat_SP,  
		 //				   para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelBGKCompSP27 execution failed");
		 //KernelKumNewCompSpongeSP27( para->getParD(0)->numberofthreads,       
			//						 para->getParD(0)->omega, 
			//						 para->getParD(0)->geoSP, 
			//						 para->getParD(0)->neighborX_SP, 
			//						 para->getParD(0)->neighborY_SP, 
			//						 para->getParD(0)->neighborZ_SP,
			//						 para->getParD(0)->coordX_SP,
			//						 para->getParD(0)->coordY_SP,
			//						 para->getParD(0)->coordZ_SP,
			//						 para->getParD(0)->d0SP.f[0],    
			//						 para->getParD(0)->size_Mat_SP,
			//						 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCasSPKum27 execution failed");
		 //printf("Level: %d \n", 0);
		 //KernelKumNewCompSP27(   para->getParD(0)->numberofthreads,       
			//					 para->getParD(0)->omega,			// omegaSponge, //
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->size_Mat_SP,
			//					 para->getParD(0)->size_Array_SP,
			//					 0,
			//					 para->getForcesDev(),
			//					 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCasSPKum27 execution failed");
		 //////////////////////////////////////////////////////////////////////////
		 //KernelKumCompSP27(  para->getParD(0)->numberofthreads,       
			//				 para->getParD(0)->omega, 
			//				 para->getParD(0)->geoSP, 
			//				 para->getParD(0)->neighborX_SP, 
			//				 para->getParD(0)->neighborY_SP, 
			//				 para->getParD(0)->neighborZ_SP,
			//				 para->getParD(0)->d0SP.f[0],    
			//				 para->getParD(0)->size_Mat_SP,  
			//				 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCasSPKum27 execution failed");
		 //////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////////
		//Wale 
		if (para->getUseWale())
		{
			KernelWaleCumOneCompSP27(para->getParD(0)->numberofthreads,
									 para->getParD(0)->omega,			
									 para->getParD(0)->geoSP, 
									 para->getParD(0)->neighborX_SP, 
									 para->getParD(0)->neighborY_SP, 
									 para->getParD(0)->neighborZ_SP,
									 para->getParD(0)->neighborWSB_SP,
								     para->getParD(0)->vx_SP,        
								     para->getParD(0)->vy_SP,        
								     para->getParD(0)->vz_SP,        
									 para->getParD(0)->d0SP.f[0],
									 para->getParD(0)->turbViscosity,
									 para->getParD(0)->size_Mat_SP,
									 para->getParD(0)->size_Array_SP,
									 0,
									 para->getForcesDev(),
									 para->getParD(0)->evenOrOdd); 
			getLastCudaError("KernelWaleCumOneCompSP27 execution failed");

			//KernelWaleCumAA2016CompSP27(para->getParD(0)->numberofthreads,
			//							para->getParD(0)->omega,			
			//							para->getParD(0)->geoSP, 
			//							para->getParD(0)->neighborX_SP, 
			//							para->getParD(0)->neighborY_SP, 
			//							para->getParD(0)->neighborZ_SP,
			//							para->getParD(0)->neighborWSB_SP,
			//							para->getParD(0)->vx_SP,        
			//							para->getParD(0)->vy_SP,        
			//							para->getParD(0)->vz_SP,        
			//							para->getParD(0)->d0SP.f[0],
			//							para->getParD(0)->turbViscosity,
			//							para->getParD(0)->size_Mat_SP,
			//							para->getParD(0)->size_Array_SP,
			//							0,
			//							para->getForcesDev(),
			//							para->getParD(0)->evenOrOdd); 
			//getLastCudaError("KernelWaleCumAA2016CompSP27 execution failed");

		} 
		else
		{
			KernelKumNewCompSP27(para->getParD(0)->numberofthreads,       
								 para->getParD(0)->omega,			
								 para->getParD(0)->geoSP, 
								 para->getParD(0)->neighborX_SP, 
								 para->getParD(0)->neighborY_SP, 
								 para->getParD(0)->neighborZ_SP,
								 para->getParD(0)->d0SP.f[0],    
								 para->getParD(0)->size_Mat_SP,
								 para->getParD(0)->size_Array_SP,
								 0,
								 para->getForcesDev(),
								 para->getParD(0)->evenOrOdd); 
			getLastCudaError("KernelCasSPKum27 execution failed");

			//KernelKumAA2016CompSP27(para->getParD(0)->numberofthreads,       
			//					 para->getParD(0)->omega, 
			//					 para->getParD(0)->geoSP, 
			//					 para->getParD(0)->neighborX_SP, 
			//					 para->getParD(0)->neighborY_SP, 
			//					 para->getParD(0)->neighborZ_SP,
			//					 para->getParD(0)->d0SP.f[0],    
			//					 para->getParD(0)->size_Mat_SP,
			//					 0,
			//					 para->getForcesDev(),
			//					 para->getParD(0)->evenOrOdd); 
			//getLastCudaError("KernelKumAA2016CompSP27 execution failed");

			//KernelCumulantD3Q27All4(para->getParD(0)->numberofthreads,
			//						 para->getParD(0)->omega, 
			//						 para->getParD(0)->geoSP, 
			//						 para->getParD(0)->neighborX_SP, 
			//						 para->getParD(0)->neighborY_SP, 
			//						 para->getParD(0)->neighborZ_SP,
			//						 para->getParD(0)->d0SP.f[0],    
			//						 para->getParD(0)->size_Mat_SP,
			//						 0,
			//						 para->getForcesDev(),
			//						 para->getParD(0)->evenOrOdd); 
			//getLastCudaError("KernelCumulantD3Q27All4 execution failed");
		}
		//////////////////////////////////////////////////////////////////////////


		////////////////////////////////////////////////////////////////////////////
		////porous media
		//if (para->getSimulatePorousMedia())
		//{
		//	KernelPMCumOneCompSP27( para->getParD(0)->numberofthreads,
		//							para->getParD(0)->omega,			
		//							para->getParD(0)->neighborX_SP, 
		//							para->getParD(0)->neighborY_SP, 
		//							para->getParD(0)->neighborZ_SP,
		//							para->getParD(0)->d0SP.f[0],    
		//							para->getParD(0)->size_Mat_SP,
		//							0,
		//							para->getForcesDev(),
		//							pm0->getPorosity(),
		//							pm0->getDarcyLBM(),
		//							pm0->getForchheimerLBM(),
		//							pm0->getSizePM(),
		//							pm0->getHostNodeIDsPM(),
		//							para->getParD(0)->evenOrOdd); 
		//	getLastCudaError("KernelPMCumOneCompSP27 execution failed");
		//}
		////////////////////////////////////////////////////////////////////////////




		 //////////////////////////////////////////////////////////////////////////
		 //incomp
 		 //KernelBGKPlusSP27(para->getParD(0)->numberofthreads,       
 		 //				   para->getParD(0)->omega, 
 		 //				   para->getParD(0)->geoSP, 
 		 //				   para->getParD(0)->neighborX_SP, 
 		 //				   para->getParD(0)->neighborY_SP, 
 		 //				   para->getParD(0)->neighborZ_SP,
 		 //				   para->getParD(0)->d0SP.f[0],    
 		 //				   para->getParD(0)->size_Mat_SP,  
 		 //				   para->getParD(0)->evenOrOdd); 
 		 //getLastCudaError("KernelBGKPlusSP27 execution failed");
		 //KernelBGKSP27(para->getParD(0)->numberofthreads,       
		 //			   para->getParD(0)->omega, 
		 //			   para->getParD(0)->geoSP, 
		 //			   para->getParD(0)->neighborX_SP, 
		 //			   para->getParD(0)->neighborY_SP, 
		 //			   para->getParD(0)->neighborZ_SP,
		 //			   para->getParD(0)->d0SP.f[0],    
		 //			   para->getParD(0)->size_Mat_SP,  
		 //			   para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelBGKP27 execution failed");
//  		 KernelMRTSP27(   para->getParD(0)->numberofthreads,       
//  							 para->getParD(0)->omega, 
//  							 para->getParD(0)->geoSP, 
//  							 para->getParD(0)->neighborX_SP, 
//  							 para->getParD(0)->neighborY_SP, 
//  							 para->getParD(0)->neighborZ_SP,
//  							 para->getParD(0)->d0SP.f[0],    
//  							 para->getParD(0)->size_Mat_SP,  
//  							 para->getParD(0)->evenOrOdd); 
//  		 getLastCudaError("KernelMRT27 execution failed");
// 		KernelCascadeSP27(  para->getParD(0)->numberofthreads,       
// 							 para->getParD(0)->omega, 
// 							 para->getParD(0)->geoSP, 
// 							 para->getParD(0)->neighborX_SP, 
// 							 para->getParD(0)->neighborY_SP, 
// 							 para->getParD(0)->neighborZ_SP,
// 							 para->getParD(0)->d0SP.f[0],    
// 							 para->getParD(0)->size_Mat_SP,  
// 							 para->getParD(0)->evenOrOdd); 
// 		 getLastCudaError("KernelCas27 execution failed");
		 //KernelKumNewSP27(   para->getParD(0)->numberofthreads,       
			//				 para->getParD(0)->omega, 
			//				 para->getParD(0)->geoSP, 
			//				 para->getParD(0)->neighborX_SP, 
			//				 para->getParD(0)->neighborY_SP, 
			//				 para->getParD(0)->neighborZ_SP,
			//				 para->getParD(0)->d0SP.f[0],    
			//				 para->getParD(0)->size_Mat_SP,  
			//				 para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCasSPKum27 execution failed");
		 //KernelCasKumSP27(para->getParD(0)->numberofthreads,       
			//              para->getParD(0)->omega, 
			//              para->getParD(0)->geoSP, 
			//              para->getParD(0)->neighborX_SP, 
			//              para->getParD(0)->neighborY_SP, 
			//              para->getParD(0)->neighborZ_SP,
			//              para->getParD(0)->d0SP.f[0],    
			//              para->getParD(0)->size_Mat_SP,  
			//              para->getParD(0)->evenOrOdd); 
		 //getLastCudaError("KernelCasSPKum27 execution failed");
         //KernelCasSPMSOHM27(  para->getParD(0)->numberofthreads,       
         //                     para->getParD(0)->omega, 
         //                     para->getParD(0)->geoSP, 
         //                     para->getParD(0)->neighborX_SP, 
         //                     para->getParD(0)->neighborY_SP, 
         //                     para->getParD(0)->neighborZ_SP,
         //                     para->getParD(0)->d0SP.f[0],    
         //                     para->getParD(0)->size_Mat_SP,  
         //                     para->getParD(0)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
         //KernelCasSPMS27(  para->getParD(0)->numberofthreads,       
         //                  para->getParD(0)->omega, 
         //                  para->getParD(0)->geoSP, 
         //                  para->getParD(0)->neighborX_SP, 
         //                  para->getParD(0)->neighborY_SP, 
         //                  para->getParD(0)->neighborZ_SP,
         //                  para->getParD(0)->d0SP.f[0],    
         //                  para->getParD(0)->size_Mat_SP,  
         //                  para->getParD(0)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
         //KernelCasSP27(para->getParD(0)->numberofthreads,       
         //              para->getParD(0)->omega, 
         //              para->getParD(0)->geoSP, 
         //              para->getParD(0)->neighborX_SP, 
         //              para->getParD(0)->neighborY_SP, 
         //              para->getParD(0)->neighborZ_SP,
         //              para->getParD(0)->d0SP.f[0],    
         //              para->getParD(0)->size_Mat_SP,  
         //              para->getParD(0)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
         //KernelCas27(para->getParD(0)->gridNX,  para->getParD(0)->gridNY,   para->getParD(0)->gridNZ, ic.s9, para->getParD(0)->geo, 
         //            para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
         //            para->getParD(0)->d0.f[0], para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd); 
         //getLastCudaError("Kernel execution failed");
         ////////////////////////////////////////////////////////////////////////////////
		 //output << "fertig \n";
		 ////////////////////////////////////////////////////////////////////////////////
         if (para->getDiffOn()==true)
         {
            if (para->getDiffMod() == 7)
            {
				//output << " Diff Mod 7\n";
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // incomp
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelADincomp7(para->getParD(0)->numberofthreads,    para->getParD(0)->diffusivity,  para->getParD(0)->geoSP, 
							   para->getParD(0)->neighborX_SP,       para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
							   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],      para->getParD(0)->size_Mat_SP,  
							   para->getParD(0)->evenOrOdd); 
			   getLastCudaError("KernelADincomp7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion boundary condition
               QNoSlipADincompDev7( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
									para->getParD(0)->d0SP.f[0],             para->getParD(0)->d7.f[0],      para->getParD(0)->Temp.temp,  
									para->getParD(0)->diffusivity,           para->getParD(0)->Temp.k,       para->getParD(0)->QGeom.q27[0], 
									para->getParD(0)->Temp.kTemp,            para->getParD(0)->Temp.kTemp,   para->getParD(0)->omega,
									para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
									para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
               getLastCudaError("QNoSlipADincompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
			   if (t<15580)//(t>500000 && t<515580)//(t>300000 && t<315580)
			   {
                 QADVeloIncompDev7(para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
								   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],			para->getParD(0)->TempVel.tempPulse, 
								   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
								   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,	para->getParD(0)->TempVel.kTemp,  
								   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
								   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
				 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
			   else
			   {
                 QADVeloIncompDev7(para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
								   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],			para->getParD(0)->TempVel.temp, 
								   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
								   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,	para->getParD(0)->TempVel.kTemp,  
								   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
								   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
                 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
           //    QADVeloIncompDev7(  para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
								   //para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],			para->getParD(0)->TempVel.tempPulse, 
								   //para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
								   //para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,	para->getParD(0)->TempVel.kTemp,  
								   //para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
								   //para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
           //    getLastCudaError("QADVeloIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
               QADPressIncompDev7(   para->getParD(0)->numberofthreads,  para->getParD(0)->nx,				para->getParD(0)->ny,
									 para->getParD(0)->d0SP.f[0],        para->getParD(0)->d7.f[0],			para->getParD(0)->TempPress.temp, 
									 para->getParD(0)->TempPress.velo,   para->getParD(0)->diffusivity,		para->getParD(0)->TempPress.k,
									 para->getParD(0)->QPress.q27[0],    para->getParD(0)->TempPress.kTemp, para->getParD(0)->TempPress.kTemp,  
									 para->getParD(0)->omega,            para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
									 para->getParD(0)->neighborZ_SP,     para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
               getLastCudaError("QADPressIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			   

			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// comp
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion kernel
      //         KernelThS7( para->getParD(0)->numberofthreads,    para->getParD(0)->diffusivity,  para->getParD(0)->geoSP, 
      //                     para->getParD(0)->neighborX_SP,       para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //                     para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],      para->getParD(0)->size_Mat_SP,  
      //                     para->getParD(0)->evenOrOdd); 
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         QADDev7( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
      //                  para->getParD(0)->d0SP.f[0],             para->getParD(0)->d7.f[0],      para->getParD(0)->Temp.temp,  
      //                  para->getParD(0)->diffusivity,           para->getParD(0)->Temp.k,       para->getParD(0)->QGeom.q27[0], 
      //                  para->getParD(0)->Temp.kTemp,            para->getParD(0)->Temp.kTemp,   para->getParD(0)->omega,
      //                  para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //                  para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
      //         getLastCudaError("QADDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADVelDev7( para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
      //                     para->getParD(0)->d0SP.f[0],          para->getParD(0)->d7.f[0],			para->getParD(0)->TempVel.temp, 
      //                     para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
      //                     para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,     para->getParD(0)->TempVel.kTemp,  
      //                     para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,		para->getParD(0)->neighborY_SP, 
      //                     para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
      //         getLastCudaError("QADVelDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADPressDev7( para->getParD(0)->numberofthreads,  para->getParD(0)->nx,				para->getParD(0)->ny,
      //                       para->getParD(0)->d0SP.f[0],        para->getParD(0)->d7.f[0],			para->getParD(0)->TempPress.temp, 
      //                       para->getParD(0)->TempPress.velo,   para->getParD(0)->diffusivity,		para->getParD(0)->TempPress.k,
      //                       para->getParD(0)->QPress.q27[0],    para->getParD(0)->TempPress.kTemp,   para->getParD(0)->TempPress.kTemp,  
      //                       para->getParD(0)->omega,            para->getParD(0)->neighborX_SP,		para->getParD(0)->neighborY_SP, 
      //                       para->getParD(0)->neighborZ_SP,     para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
      //         getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            } 
            else if (para->getDiffMod() == 27)
            {
			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// incomp
			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelADincomp27(   para->getParD(0)->numberofthreads,    para->getParD(0)->diffusivity,  para->getParD(0)->geoSP, 
								   para->getParD(0)->neighborX_SP,       para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
								   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],     para->getParD(0)->size_Mat_SP,  
								   para->getParD(0)->evenOrOdd); 
			   getLastCudaError("KernelADincomp27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
          //     //advection diffusion boundary condition
          //     QNoSlipADincompDev27( para->getParD(0)->numberofthreads,      para->getParD(0)->nx,           para->getParD(0)->ny,
									 //para->getParD(0)->d0SP.f[0],            para->getParD(0)->d27.f[0],     para->getParD(0)->Temp.temp,  
									 //para->getParD(0)->diffusivity,          para->getParD(0)->Temp.k,       para->getParD(0)->QGeom.q27[0], 
									 //para->getParD(0)->Temp.kTemp,           para->getParD(0)->Temp.kTemp,   para->getParD(0)->omega,
									 //para->getParD(0)->neighborX_SP,         para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
									 //para->getParD(0)->size_Mat_SP,          para->getParD(0)->evenOrOdd);
          //     getLastCudaError("QNoSlipADincompDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
			   if (t>1600000 && t<1662320)//(t>400000 && t<415580)//(t>500000 && t<515580)//(t>300000 && t<315580)//(t<15580)//
			   {
				   QADVeloIncompDev27( para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
									   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],		para->getParD(0)->TempVel.tempPulse, 
									   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
									   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,   para->getParD(0)->TempVel.kTemp,  
									   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
									   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
				   getLastCudaError("QADVeloIncompDev27 execution failed");
			   }
			   else
			   {
				   QADVeloIncompDev27( para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
									   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],		para->getParD(0)->TempVel.temp, 
									   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
									   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,   para->getParD(0)->TempVel.kTemp,  
									   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
									   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
				   getLastCudaError("QADVeloIncompDev27 execution failed");
			   }
           //    QADVeloIncompDev27( para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
								   //para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],		para->getParD(0)->TempVel.temp, 
								   //para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->TempVel.k,
								   //para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->TempVel.kTemp,   para->getParD(0)->TempVel.kTemp,  
								   //para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
								   //para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
           //    getLastCudaError("QADVeloIncompDev27 execution failed");
               //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + pressure boundary condition
           //    QADPressIncompDev27(   para->getParD(0)->numberofthreads,  para->getParD(0)->nx,					para->getParD(0)->ny,
									  //para->getParD(0)->d0SP.f[0],        para->getParD(0)->d27.f[0],			para->getParD(0)->TempPress.temp, 
									  //para->getParD(0)->TempPress.velo,   para->getParD(0)->diffusivity,		para->getParD(0)->TempPress.k,
									  //para->getParD(0)->QPress.q27[0],    para->getParD(0)->TempPress.kTemp,    para->getParD(0)->TempPress.kTemp,  
									  //para->getParD(0)->omega,            para->getParD(0)->neighborX_SP,		para->getParD(0)->neighborY_SP, 
									  //para->getParD(0)->neighborZ_SP,     para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
           //    getLastCudaError("QADPressIncompDev27 execution failed");
               //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //////// comp
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion kernel
      //         KernelThS27(para->getParD(0)->numberofthreads,
      //                     para->getParD(0)->diffusivity, 
      //                     para->getParD(0)->geoSP, 
      //                     para->getParD(0)->neighborX_SP, 
      //                     para->getParD(0)->neighborY_SP, 
      //                     para->getParD(0)->neighborZ_SP,
      //                     para->getParD(0)->d0SP.f[0],    
      //                     para->getParD(0)->d27.f[0],    
      //                     para->getParD(0)->size_Mat_SP,  
      //                     para->getParD(0)->evenOrOdd); 
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         //QADBBDev27(para->getParD(0)->numberofthreads,      para->getParD(0)->nx,           para->getParD(0)->ny,
      //         //           para->getParD(0)->d0SP.f[0],            para->getParD(0)->d27.f[0],     para->getParD(0)->Temp.temp,  
      //         //           para->getParD(0)->diffusivity,          para->getParD(0)->Temp.k,       para->getParD(0)->QGeom.q27[0], 
      //         //           para->getParD(0)->Temp.kTemp,           para->getParD(0)->Temp.kTemp,   para->getParD(0)->omega,
      //         //           para->getParD(0)->neighborX_SP,         para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //         //           para->getParD(0)->size_Mat_SP,          para->getParD(0)->evenOrOdd);
      //         //getLastCudaError("QADBBDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
			   //if (t>400000 && t<415580)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
			   //{
				  // QADVelDev27(para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
						//	   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],		para->getParD(0)->TempVel.tempPulse, 
						//	   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->Qinflow.k,
						//	   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,  
						//	   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
						//	   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
			   //else
			   //{
				  // QADVelDev27(para->getParD(0)->numberofthreads,    para->getParD(0)->nx,				para->getParD(0)->ny,
						//	   para->getParD(0)->d0SP.f[0],          para->getParD(0)->d27.f[0],		para->getParD(0)->TempVel.temp, 
						//	   para->getParD(0)->TempVel.velo,       para->getParD(0)->diffusivity,		para->getParD(0)->Qinflow.k,
						//	   para->getParD(0)->Qinflow.q27[0],     para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,  
						//	   para->getParD(0)->omega,              para->getParD(0)->neighborX_SP,	para->getParD(0)->neighborY_SP, 
						//	   para->getParD(0)->neighborZ_SP,       para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         //QADPressDev27( para->getParD(0)->numberofthreads,  para->getParD(0)->nx,				para->getParD(0)->ny,
      //         //               para->getParD(0)->d0SP.f[0],        para->getParD(0)->d27.f[0],			para->getParD(0)->TempPress.temp, 
      //         //               para->getParD(0)->TempPress.velo,   para->getParD(0)->diffusivity,		para->getParD(0)->TempPress.k,
      //         //               para->getParD(0)->QPress.q27[0],    para->getParD(0)->TempPress.kTemp,  para->getParD(0)->TempPress.kTemp,  
      //         //               para->getParD(0)->omega,            para->getParD(0)->neighborX_SP,		para->getParD(0)->neighborY_SP, 
      //         //               para->getParD(0)->neighborZ_SP,     para->getParD(0)->size_Mat_SP,		para->getParD(0)->evenOrOdd);
      //         //getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
         }
		//////////////////////////////////////////////////////////////////////////////
		if (para->getNumprocs() > 1)
		{
			////1D domain decomposition
			//exchangePostCollDataGPU27(para, comm, 0);
			//////////////////////////////////////////////////////////////////////////
			//3D domain decomposition
			//output << "start exchange Post X (level 0) \n";
			exchangePostCollDataXGPU27(para.get(), comm.get(), 0);
			//output << "end exchange Post X (level 0) \n";
			//output << "start exchange Post Y (level 0) \n";
			exchangePostCollDataYGPU27(para.get(), comm.get(), 0);
			//output << "end exchange Post Y (level 0) \n";
			//output << "start exchange Post Z (level 0) \n";
			exchangePostCollDataZGPU27(para.get(), comm.get(), 0);
			//output << "end exchange Post Z (level 0) \n";
			////////////////////////////////////////////////////////////////////////
			//3D domain decomposition convection diffusion
			if (para->getDiffOn()==true)
			{
				exchangePostCollDataADXGPU27(para.get(), comm.get(), 0);
				exchangePostCollDataADYGPU27(para.get(), comm.get(), 0);
				exchangePostCollDataADZGPU27(para.get(), comm.get(), 0);
			}
		}
		////////////////////////////////////////////////////////////////////////////////
		//PressX0
		//QPressDevOld27(	para->getParD(0)->numberofthreads, para->getParD(0)->QpressX0.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QpressX0.k,  
		//				para->getParD(0)->QpressX0.kN,     para->getParD(0)->QpressX0.kQ,  para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
		//////////////////////////////////////////////////////////////////////////////
		//PressX1
		//QPressDevOld27(	para->getParD(0)->numberofthreads, para->getParD(0)->QpressX1.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QpressX1.k,  
		//				para->getParD(0)->QpressX1.kN,     para->getParD(0)->QpressX1.kQ,  para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
        ////////////////////////////////////////////////////////////////////////////////
		//output << "vor der DruckRB\n";
		//Ship
		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		//				para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
        ////////////////////////////////////////////////////////////////////////////////
		//Andrea - Soeren
		//QPressDevZero27(para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],       
		//				para->getParD(0)->QPress.k,  	   para->getParD(0)->QPress.kQ,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		//				para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
        ////////////////////////////////////////////////////////////////////////////////
			//QDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			//		para->getParD(0)->d0SP.f[0],             para->getParD(0)->QWall.k,		 para->getParD(0)->QWall.q27[0], 
			//		para->getParD(0)->kQ,                    para->getParD(0)->kQ,           para->getParD(0)->omega,
			//		para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//		para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			//getLastCudaError("QDev27 execution failed");
		//???????????????????????????????????????????
		//WallFuncDev27(  para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
		//				para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
		//				para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("WallFuncDev27 execution failed");
		//WallFuncDev27(  para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
		//				para->getParD(0)->QGeom.Vx,        para->getParD(0)->QGeom.Vy,     para->getParD(0)->QGeom.Vz,
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
		//				para->getParD(0)->QGeom.kQ,        para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("WallFuncDev27 execution failed");
		//???????????????????????????????????????????
		////////////////////////////////////////////////////////////////////////////////
		      //QVelDevice1h27(   para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
								//para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
								//para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
								//para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,          
								//para->getPhi(),                    para->getAngularVelocity(),
								//para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
								//para->getParD(0)->coordX_SP,       para->getParD(0)->coordY_SP,    para->getParD(0)->coordZ_SP,
								//para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		      //getLastCudaError("QVelDev27 execution failed");
		//???????????????????????????????????????????
		    //  QVelDev27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
		    //            para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
		    //            para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
						//para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
		    //            para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		    //            para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		    //  getLastCudaError("QVelDev27 execution failed");
		     // QVelDevComp27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
							//para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
							//para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
							//para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
							//para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
							//para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		     // getLastCudaError("QVelDevComp27 execution failed");
		     // QVelDevComp27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
							//para->getParD(0)->QGeom.Vx,        para->getParD(0)->QGeom.Vy,     para->getParD(0)->QGeom.Vz,
							//para->getParD(0)->d0SP.f[0],       para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
							//para->getParD(0)->QGeom.kQ,        para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
							//para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
							//para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		     // getLastCudaError("QVelDevComp27 execution failed");
		      //QVelDevCompZeroPress27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,             para->getParD(0)->ny,
								//	 para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,     para->getParD(0)->Qinflow.Vz,
								//	 para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,      para->getParD(0)->Qinflow.q27[0], 
								//	 para->getParD(0)->kInflowQ,        para->getParD(0)->Qinflow.kArray, para->getParD(0)->omega,
								//	 para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP,   para->getParD(0)->neighborZ_SP,
								//	 para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		      //getLastCudaError("QVelDevComp27 execution failed");

			  ////////////////////////////////////////////////////////////////////////////
		      //QVeloDevEQ27(para->getParD(0)->numberofthreads,
						  // para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
						  // para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k, 
						  // para->getParD(0)->kInflowQ,        para->getParD(0)->omega,
						  // para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
						  // para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		      //getLastCudaError("QVeloDevEQ27 execution failed");
			  ////////////////////////////////////////////////////////////////////////////


			//  ////////////////////////////////////////////////////////////////////////////
		 //   if ((para->getParD(0)->QGeom.kQ) > 0 && (para->getIsGeometryValues()))
			//{
			//  //QVelDevCompZeroPress27(	para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
			//		//					para->getParD(0)->QGeom.Vx,        para->getParD(0)->QGeom.Vy,     para->getParD(0)->QGeom.Vz,
			//		//					para->getParD(0)->d0SP.f[0],       para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
			//  //							para->getParD(0)->QGeom.kQ,        para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
			//		//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//		//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//  //getLastCudaError("QVelDevCompZeroPress27 execution failed");
			//  //////////////////////////////////////////////////////////////////////////
			//   QVelDevCompPlusSlip27( para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
			//						  para->getParD(0)->QGeom.Vx,        para->getParD(0)->QGeom.Vy,     para->getParD(0)->QGeom.Vz,
			//						  para->getParD(0)->d0SP.f[0],       para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
			//						  para->getParD(0)->QGeom.kQ,        para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
			//						  para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//						  para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//   getLastCudaError("QVelDevCompPlusSlip27 execution failed");

			//}

		////////////////////////////////////////////////////////////////////////////////
		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		//				para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QOutletPressDev27 execution failed");
	    ////////////////////////////////////////////////////////////////////////////////
		//QPressDevZero27(para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],       
		//				para->getParD(0)->QPress.k,  	   para->getParD(0)->QPress.kQ,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
        ////////////////////////////////////////////////////////////////////////////////
		//Inlet - Outlet
		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QInlet.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QInlet.k,  
		//				para->getParD(0)->QInlet.kN,       para->getParD(0)->QInlet.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QInletPressDev27 execution failed");
		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QOutlet.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QOutlet.k,  
		//				para->getParD(0)->QOutlet.kN,      para->getParD(0)->QOutlet.kQ,   para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QOutletPressDev27 execution failed");
		//QPressDevZero27(para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],       
		//				para->getParD(0)->QOutlet.k,  	   para->getParD(0)->QOutlet.kQ,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");
		////////////////////////////////////////////////////////////////////////////////
         //if (  myid == numprocs - 1)  
		//QPressDev27_IntBB(  para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC,
		//					para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,       para->getParD(0)->QPress.q27[0], 
		//					para->getParD(0)->QPress.kQ,       para->getParD(0)->QPress.kQ,      para->getParD(0)->omega,
		//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP,   para->getParD(0)->neighborZ_SP,
		//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27_IntBB X1 coarse execution failed");
         //   QPressDevOld27(para->getParD(0)->numberofthreads, para->getParD(0)->Qoutflow.RhoBC, 
         //                  para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qoutflow.k,  
						   //para->getParD(0)->Qoutflow.kN,     para->getParD(0)->kOutflowQ,    para->getParD(0)->omega,
         //                  para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
         //                  para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
         //   getLastCudaError("QPressDev27 execution failed");
            //QPressDevDirDepBot27( para->getParD(0)->numberofthreads,       RhoBCOutflowD,
            //                      para->getParD(0)->d0SP.f[0],    QoutflowD.k, kOutflowQ, ic.s9,
            //                      para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
            //                      para->getParD(0)->size_Mat_SP,  para->getParD(0)->evenOrOdd);
            //getLastCudaError("QPressDev27 execution failed");
            //QPressDevFixBackflow27( para->getParD(0)->numberofthreads,       RhoBCOutflowD,
            //                        para->getParD(0)->d0SP.f[0],    QoutflowD.k, kOutflowQ, ic.s9,
            //                        para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
            //                        para->getParD(0)->size_Mat_SP,  para->getParD(0)->evenOrOdd);
            //getLastCudaError("QPressDev27 execution failed");
            //   QPressDev27(para->getParD(0)->numberofthreads,   para->getParD(0)->nx, para->getParD(0)->ny,
            //               RhoBCOutflowD,
            //               para->getParD(0)->d0SP.f[0],  QoutflowD.k, QoutflowD.q27[0], kOutflowQ, kOutflowQ, ic.s9,
            //               para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
            //               para->getParD(0)->size_Mat_SP, para->getParD(0)->evenOrOdd);
            //getLastCudaError("QPressDev27 execution failed");
            //   QPressDev27(para->getParD(0)->numberofthreads,   para->getParD(0)->nx, para->getParD(0)->ny,
            //               RhoBCOutflowD,
            //               para->getParD(0)->d0.f[0],  QoutflowD.k, QoutflowD.q27[0], kOutflowQ, kOutflowQ, ic.s9,
            //               para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
            //               para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd);
            //getLastCudaError("Kernel execution failed");
            //   BcPress27( para->getParD(0)->nx,  para->getParD(0)->ny,      para->getParD(0)->gridNZ-1, para->getParD(0)->gridNX,  para->getParD(0)->gridNY, 
            //   para->getParD(0)->geo,     para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
            //   para->getParD(0)->d0.f[0], para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd );
            //getLastCudaError("Kernel execution failed");
         ////////////////////////////////////////////////////////////////////////////////

         ////////////////////////////////////////////////////////////////////////////////
         //if (  myid == 0) 
         //   QVelDeviceCouhette27(para->getParD(0)->numberofthreads, para->getParD(0)->Qinflow.Vx,   para->getParD(0)->Qinflow.Vy,     para->getParD(0)->Qinflow.Vz,
								 //para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
								 //para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
								 //para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
								 //para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
         //   getLastCudaError("QVelDeviceCouhette27 execution failed");
        //    QVelDevicePlainBB27(para->getParD(0)->numberofthreads, para->getParD(0)->Qinflow.Vx,   para->getParD(0)->Qinflow.Vy,     para->getParD(0)->Qinflow.Vz,
								//para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
								//para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
								//para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
								//para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
        //    getLastCudaError("QVelDevPlainBB27 execution failed");
      //      QVelDev27(  para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
      //                  para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
      //                  para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
						//para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
      //                  para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //                  para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
      //      getLastCudaError("QVelDev27 execution failed");
            //QVelDev27(para->getParD(0)->numberofthreads,   para->getParD(0)->nx, para->getParD(0)->ny,
            //          VxInflowD, VyInflowD, VzInflowD,
            //          para->getParD(0)->d0.f[0],  QinflowD.k, QinflowD.q27[0], kInflowQ, kInflowQ, ic.s9,
            //          para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
            //          para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd);
            //          getLastCudaError("Kernel execution failed");
            //QDev27(para->getParD(0)->numberofthreads,   para->getParD(0)->nx, para->getParD(0)->ny,
            //             para->getParD(0)->d0.f[0],  QinflowD.k, QinflowD.q27[0], kInflowQ, kInflowQ, ic.s9,
            //             para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
            //             para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd);
            //          getLastCudaError("Kernel execution failed");
            //BBDev27(para->getParD(0)->numberofthreads,   para->getParD(0)->nx, para->getParD(0)->ny,
            //             para->getParD(0)->d0.f[0],  QinflowD.k, QinflowD.q27[0], kInflowQ, kInflowQ, ic.s9,
            //             para->getParD(0)->neighborX, para->getParD(0)->neighborY, para->getParD(0)->neighborZ,
            //             para->getParD(0)->size_Mat, para->getParD(0)->evenOrOdd);
            //          getLastCudaError("Kernel execution failed");
           ////////////////////////////////////////////////////////////////////////////
		   //Slip
			//at the wall
		   //QSlipDevComp27(	para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],    para->getParD(0)->QWall.k,
		   //					para->getParD(0)->QWall.q27[0],    para->getParD(0)->kQ,		   para->getParD(0)->omega,
		   //					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		   //					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		   //getLastCudaError("QSlipDev27 execution failed");
		 //  if (para->getParD(0)->kSlipQ > 0)
			//{
			//	QSlipDevComp27( para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],    para->getParD(0)->QSlip.k,
			//					para->getParD(0)->QSlip.q27[0],    para->getParD(0)->kSlipQ,       para->getParD(0)->omega,
			//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//	getLastCudaError("QSlipDev27 execution failed");
			//	//QDevComp27(para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			//	//			 para->getParD(0)->d0SP.f[0],             para->getParD(0)->QSlip.k,      para->getParD(0)->QSlip.q27[0],
			//	//			 para->getParD(0)->QSlip.kQ,              para->getParD(0)->QSlip.kQ,     para->getParD(0)->omega,
			//	//			 para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//	//			 para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			//	//getLastCudaError("QDevComp27 (Geom) execution failed");
		 //  }
			//QSlipDev27( para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],    para->getParD(0)->QSlip.k,
			//			para->getParD(0)->QSlip.q27[0],    para->getParD(0)->kSlipQ,       para->getParD(0)->omega,
			//			para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//			para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//getLastCudaError("QSlipDev27 execution failed");
			//QSlipDev27( para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],    para->getParD(0)->QGeom.k,
			//			para->getParD(0)->QGeom.q27[0],    para->getParD(0)->QGeom.kQ,       para->getParD(0)->omega,
			//			para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//			para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//getLastCudaError("QSlipDev27 execution failed");
			//////////////////////////////////////////////////////////////////////////////
			//output << "vor der WallBC\n";
			////Test 2ndOrder
            //BBDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
            //         para->getParD(0)->d0SP.f[0],             QD.k, QD.q27[0], kQ, kQ,        para->getParD(0)->omega,
            //         para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
            //         para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
            //getLastCudaError("QDev27 execution failed");
			//////////////////////////////////////////////////////////////////////////
			//Wall
      //      BBDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
      //               para->getParD(0)->d0SP.f[0],             para->getParD(0)->QWall.k,      para->getParD(0)->QWall.q27[0], 
					 //para->getParD(0)->kQ,                    para->getParD(0)->kQ,           para->getParD(0)->omega,
      //               para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //               para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
      //      getLastCudaError("QDev27 execution failed");
			//QDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			//		para->getParD(0)->d0SP.f[0],             para->getParD(0)->QWall.k,		 para->getParD(0)->QWall.q27[0], 
			//		para->getParD(0)->kQ,                    para->getParD(0)->kQ,           para->getParD(0)->omega,
			//		para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//		para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			//getLastCudaError("QDev27 execution failed");
		   //if (para->getParD(0)->kQ > 0)
		   //{
			  // QDevComp27(para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			  // 		      para->getParD(0)->d0SP.f[0],             para->getParD(0)->QWall.k,	   para->getParD(0)->QWall.q27[0], 
			  // 		      para->getParD(0)->kQ,                    para->getParD(0)->kQ,           para->getParD(0)->omega,
			  // 		      para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			  // 		      para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			  // getLastCudaError("QDevComp27 (Wall) execution failed");
		   //}
			////////////////////////////////////////////////////////////////////////////
			//Geom
      //      BBDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
      //               para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
					 //para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
      //               para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
      //               para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
      //      getLastCudaError("QDev27 execution failed");
			//QDev27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			//		para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,		 para->getParD(0)->QGeom.q27[0], 
			//		para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
			//		para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//		para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			//getLastCudaError("QDev27 execution failed");
			//QDevComp27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
			//			para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,		 para->getParD(0)->QGeom.q27[0], 
			//			para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
			//			para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//			para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
			//getLastCudaError("QDevComp27 (Geom) execution failed");


		//QPressDevOld27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		//				para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");

		//QPressDevEQZ27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		//				para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDevEQZ27 execution failed");

		//QPressDevZero27(para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],       
		//				para->getParD(0)->QPress.k,  	   para->getParD(0)->QPress.kQ,
		//				para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QPressDev27 execution failed");

		  //QDevComp27(para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
		  //			 para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
		  //			 para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
		  //			 para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		  //			 para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QDevComp27 (Geom) execution failed");


			//////////////////////////////////////////////////////////////////////////////////
			////print forces to file
			//////printForcing(para);
			////if(para->getTOut()>0 && t%para->getTOut()==0 && t>para->getStartTurn())              printForcing(para);
			//////////////////////////////////////////////////////////////////////////////////
			////get velocities at slip bc to set the forcing 
			//GetVelotoForce27(para->getParD(0)->numberofthreads, para->getParD(0)->d0SP.f[0],    para->getParD(0)->QSlip.k, 
			//				 para->getParD(0)->kSlipQ,          
			//				 para->getParD(0)->VxForce,         para->getParD(0)->VyForce,      para->getParD(0)->VzForce,
			//	             para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//	             para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//////////////////////////////////////////////////////////////////////////////////
			////calculate the new forcing
			//calcVeloForce(para);
			//////////////////////////////////////////////////////////////////////////////////

          ////////////////////////////////////////////////////////////////////////////////
          if (para->getParD(0)->evenOrOdd==true)  para->getParD(0)->evenOrOdd=false;
          else                                    para->getParD(0)->evenOrOdd=true;
          ////////////////////////////////////////////////////////////////////////////////


		if (para->getUseWale())
		{
		    CalcMacCompSP27(para->getParD(0)->vx_SP,       
						    para->getParD(0)->vy_SP,        
						    para->getParD(0)->vz_SP,        
						    para->getParD(0)->rho_SP, 
						    para->getParD(0)->press_SP, 
						    para->getParD(0)->geoSP,       
						    para->getParD(0)->neighborX_SP, 
						    para->getParD(0)->neighborY_SP, 
						    para->getParD(0)->neighborZ_SP,
						    para->getParD(0)->size_Mat_SP,
						    para->getParD(0)->numberofthreads,       
						    para->getParD(0)->d0SP.f[0],    
						    para->getParD(0)->evenOrOdd);
            getLastCudaError("CalcMacSP27 execution failed"); 

		} 
		  //////////////////////////////////////////////////////////////////////////////////
		  ////calculate the new forcing
		  //if (((t%10) == 0) && (t >= 10)/*((t%para->getTStartOut()) == 0) && (t >= para->getTStartOut())*/)
		  //{
			 // forceCalculator->calcPIDControllerForForce(para);
		  //}
		  //////////////////////////////////////////////////////////////////////////////////
		  
		  // ////////////////////////////////////////////////////////////////////////
		 // if (para->getDiffOn()==true)
		 // {
			//  if (para->getDiffMod() == 7)
			//  {
			//	  //?????????????????
			//  }
			//  else if (para->getDiffMod() == 27)
			//  {
			//		QPressNoRhoDev27(   para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
			//							para->getParD(0)->d27.f[0],        para->getParD(0)->QPress.k,  
			//							para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
			//							para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//							para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//		getLastCudaError("QPressDev27 execution failed");
			//  }
		 // }
		 // ////////////////////////////////////////////////////////////////////////



		 //////////////////////////////////////////////////////////////////////////////////
		 ////press EQ comp
		 //if (para->getParD(0)->QPress.kQ > 0)
		 //{
			// // //////////////////////////////////////////////////////////////////////////////////
			// QPressNoRhoDev27(  para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
			// 					para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
			// 					para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
			// 					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			// 					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			// getLastCudaError("QPressNoRhoDev27 execution failed");
			// //QPressDevEQZ27(para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		 //	//				para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
			//	//			para->getParD(0)->QPress.kN,       para->getParD(0)->kDistTestRE.f[0],       
			//	//			para->getParD(0)->QPress.kQ,       para->getParD(0)->omega,
			//	//			para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		 //	//				para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			// //getLastCudaError("QPressDevEQZ27 execution failed");

 		////	QInflowScaleByPressDev27(   para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
			////							para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
			////							para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
			////							para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			////							para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			////getLastCudaError("QInflowScaleByPressDev27 execution failed");

		 //}
		 //////////////////////////////////////////////////////////////////////////////////
		 ////only for a round off error test
		 //para->cudaCopyTestREtoHost(0,para->getParH(0)->QPress.kQ);
		 //printRE(para, t);
		 //////////////////////////////////////////////////////////////////////////////////

		  //////////////////////////////////////////////////////////////////////////////////
			//QPressNoRhoDev27(   para->getParD(0)->numberofthreads, para->getParD(0)->Qoutflow.RhoBC, 
			//					para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qoutflow.k,  
			//					para->getParD(0)->Qoutflow.kN,     para->getParD(0)->Qoutflow.kQ,    para->getParD(0)->omega,
			//					para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
			//					para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
			//getLastCudaError("QPressNoRhoDev27 execution failed");
		  //////////////////////////////////////////////////////////////////////////////////
		  ////press NEQ incomp
		  //QPressDevIncompNEQ27(para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		  //					   para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		  //					   para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		  //					   para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		  //					   para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QPressDevIncompNEQ27 execution failed");
		  //////////////////////////////////////////////////////////////////////////////////
		  //press NEQ comp
		  //QPressDevNEQ27( para->getParD(0)->numberofthreads, para->getParD(0)->QPress.RhoBC, 
		  //				  para->getParD(0)->d0SP.f[0],       para->getParD(0)->QPress.k,  
		  //				  para->getParD(0)->QPress.kN,       para->getParD(0)->QPress.kQ,    para->getParD(0)->omega,
		  //				  para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		  //				  para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QPressDevNEQ27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////
          //if (  myid == numprocs - 1)  
          //   PressSchlaffer27( para->getParD(0)->numberofthreads,  para->getParD(0)->Qoutflow.RhoBC,
          //                     para->getParD(0)->d0SP.f[0],        para->getParD(0)->Qoutflow.Vx, 
							   //para->getParD(0)->Qoutflow.Vy,      para->getParD(0)->Qoutflow.Vz, 
							   //para->getParD(0)->Qoutflow.deltaVz, para->getParD(0)->Qoutflow.k,  
							   //para->getParD(0)->Qoutflow.kN,      para->getParD(0)->kOutflowQ,                      
							   //para->getParD(0)->omega,            para->getParD(0)->neighborX_SP,    
							   //para->getParD(0)->neighborY_SP,     para->getParD(0)->neighborZ_SP,
          //                     para->getParD(0)->size_Mat_SP,      para->getParD(0)->evenOrOdd);
          //getLastCudaError("PressSchlaffer27 execution failed");
          ////////////////////////////////////////////////////////////////////////////////
          //if (  myid == 0)  
       //      VelSchlaffer27(para->getParD(0)->numberofthreads, t,
       //                     para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.Vz, 
							//para->getParD(0)->Qinflow.deltaVz, para->getParD(0)->Qinflow.k,  
							//para->getParD(0)->Qinflow.kN,      para->getParD(0)->kInflowQ, 
							//para->getParD(0)->omega,           para->getParD(0)->neighborX_SP, 
							//para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
       //                     para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
       //   getLastCudaError("VelSchlaffer27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////
		  ////High noon incomp
		  //QVelDevIncompHighNu27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
				//			    para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
				//			    para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
				//			    para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
				//			    para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
				//			    para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QVelDevComp27 execution failed");

		  //QDevIncompHighNu27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
				//			  para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
				//			  para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
				//			  para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
				//			  para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QDevComp27 (Geom) execution failed");
		  //////////////////////////////////////////////////////////////////////////////////
		  ////High noon comp
		  //QVelDevCompHighNu27(para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
				//			  para->getParD(0)->Qinflow.Vx,      para->getParD(0)->Qinflow.Vy,   para->getParD(0)->Qinflow.Vz,
				//			  para->getParD(0)->d0SP.f[0],       para->getParD(0)->Qinflow.k,    para->getParD(0)->Qinflow.q27[0], 
				//			  para->getParD(0)->kInflowQ,        para->getParD(0)->kInflowQ,     para->getParD(0)->omega,
				//			  para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
				//			  para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QVelDevComp27 execution failed");

		  //QDevCompHighNu27( para->getParD(0)->numberofthreads,       para->getParD(0)->nx,           para->getParD(0)->ny,
		  //					para->getParD(0)->d0SP.f[0],             para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
		  //					para->getParD(0)->QGeom.kQ,              para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
		  //					para->getParD(0)->neighborX_SP,          para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		  //					para->getParD(0)->size_Mat_SP,           para->getParD(0)->evenOrOdd);
		  //getLastCudaError("QDevComp27 (Geom) execution failed");
		  //////////////////////////////////////////////////////////////////////////////////



         if (para->getMaxLevel()>=1)
         {
            //////////////////////////////////////////////////////////////////////////
            // Scaling CF and FC
            //////////////////////////////////////////////////////////////////////////
            //ScaleCF27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
            //            para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
            //            para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP, para->getParD(1)->neighborZ_SP,
            //            para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,  para->getParD(0)->evenOrOdd,
            //            para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
            //            para->getParD(0)->K_CF,           para->getParD(0)->omega,        para->getParD(1)->omega, 
            //            para->getParD(0)->vis,            para->getParD(0)->nx,           para->getParD(0)->ny, 
            //            para->getParD(1)->nx,             para->getParD(1)->ny,           para->getParD(0)->gridNX);
            //getLastCudaError("ScaleCF27 execution failed");

            //ScaleFC27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
            //            para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP, 
            //            para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP, para->getParD(1)->neighborZ_SP, 
            //            para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,  para->getParD(0)->evenOrOdd,
            //            para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
            //            para->getParD(0)->K_FC,           para->getParD(0)->omega,        para->getParD(1)->omega, 
            //            para->getParD(0)->vis,            para->getParD(0)->nx,           para->getParD(0)->ny, 
            //            para->getParD(1)->nx,             para->getParD(1)->ny,           para->getParD(0)->gridNX);
            //getLastCudaError("ScaleFC27 execution failed");
            //////////////////////////////////////////////////////////////////////////
            //ScaleCFEff27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
            //               para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offCF);
            //getLastCudaError("ScaleCF27 execution failed");

            //ScaleFCEff27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
            //               para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offFC);
            //getLastCudaError("ScaleFC27 execution failed");
            //////////////////////////////////////////////////////////////////////////
            //ScaleCFLast27( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
            //               para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offCF);
            //getLastCudaError("ScaleCF27 execution failed");

            //ScaleFCLast27( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
            //               para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offFC);
            //getLastCudaError("ScaleFC27 execution failed");
            //////////////////////////////////////////////////////////////////////////
            //ScaleCFpress27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
            //               para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offCF);
            //getLastCudaError("ScaleCF27 execution failed");

            //ScaleFCpress27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
            //               para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offFC);
            //getLastCudaError("ScaleFC27 execution failed");


            //////////////////////////////////////////////////////////////////////////
			//fine to coarse interpolation
            //ScaleFC_0817_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
							     //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
							     //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
							     //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							     //para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
							     //para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							     //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							     //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							     //para->getParD(0)->offFC);
            //getLastCudaError("ScaleFC_0817_comp_27 execution failed");
            ////////////////////////////////////////////////////////////////////////////
           // ScaleFC_Fix_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
							    //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
							    //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
							    //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							    //para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
							    //para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							    //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							    //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							    //para->getParD(0)->offFC);
           // getLastCudaError("ScaleFC_Fix_comp_27 execution failed");
            ////////////////////////////////////////////////////////////////////////////

            ScaleFC_RhoSq_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
							      para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
							      para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
							      para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							      para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
							      para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							      para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							      para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							      para->getParD(0)->offFC);
            getLastCudaError("ScaleFC_RhoSq_comp_27 execution failed");

		    //////////////////////////////////////////////////////////////////////////
           // ScaleFC_RhoSq_3rdMom_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
										 //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
										 //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
										 //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
										 //para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
										 //para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
										 //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
										 //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
										 //para->getParD(0)->offFC);
           // getLastCudaError("ScaleFC_RhoSq_3rdMom_comp_27 execution failed");
			//////////////////////////////////////////////////////////////////////////
         //   ScaleFC_AA2016_comp_27( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
									//para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
									//para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
									//para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
									//para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
									//para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
									//para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
									//para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
									//para->getParD(0)->offFC);
         //   getLastCudaError("ScaleFC_AA2016_comp_27 execution failed");
			//////////////////////////////////////////////////////////////////////////
			//data exchange
			 if (para->getNumprocs() > 1)
			 {
				 ////1D domain decomposition
				 //exchangePreCollDataGPU27(para, comm.get(), 0);
				 //3D domain decomposition
				 //output << "start exchange Pre X (level 0) \n";
				 exchangePreCollDataXGPU27(para.get(), comm.get(), 0);
				 //output << "end exchange Pre X (level 0) \n";
				 //output << "start exchange Pre Y (level 0) \n";
				 exchangePreCollDataYGPU27(para.get(), comm.get(), 0);
				 //output << "end exchange Pre Y (level 0) \n";
				 //output << "start exchange Pre Z (level 0) \n";
				 exchangePreCollDataZGPU27(para.get(), comm.get(), 0);
				 //output << "end exchange Pre Z (level 0) \n";
				 //////////////////////////////////////////////////////////////////////////
				 //3D domain decomposition convection diffusion
				 if (para->getDiffOn()==true)
				 {
					 exchangePreCollDataADXGPU27(para.get(), comm.get(), 0);
					 exchangePreCollDataADYGPU27(para.get(), comm.get(), 0);
					 exchangePreCollDataADZGPU27(para.get(), comm.get(), 0);
				 }
			 }
		    //////////////////////////////////////////////////////////////////////////
		    //coarse to fine interpolation
             //ScaleCF_0817_comp_27(para->getParD(0)->d0SP.f[0], para->getParD(1)->d0SP.f[0],
             //    para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
             //    para->getParD(1)->neighborX_SP, para->getParD(1)->neighborY_SP, para->getParD(1)->neighborZ_SP,
             //    para->getParD(0)->size_Mat_SP, para->getParD(1)->size_Mat_SP, para->getParD(0)->evenOrOdd,
             //    para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF,
             //    para->getParD(0)->K_CF, para->getParD(0)->omega, para->getParD(1)->omega,
             //    para->getParD(0)->vis, para->getParD(0)->nx, para->getParD(0)->ny,
             //    para->getParD(1)->nx, para->getParD(1)->ny, para->getParD(0)->numberofthreads,
             //    para->getParD(0)->offCF);
             //getLastCudaError("ScaleCF_0817_comp_27 execution failed");
		    ////////////////////////////////////////////////////////////////////////
           // ScaleCF_Fix_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
							    //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
							    //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
							    //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							    //para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
							    //para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							    //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							    //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							    //para->getParD(0)->offCF);
           // getLastCudaError("ScaleCF_Fix_comp_27 execution failed");
		    ////////////////////////////////////////////////////////////////////////


             ScaleCF_RhoSq_comp_27(para->getParD(0)->d0SP.f[0], para->getParD(1)->d0SP.f[0],
                 para->getParD(0)->neighborX_SP, para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
                 para->getParD(1)->neighborX_SP, para->getParD(1)->neighborY_SP, para->getParD(1)->neighborZ_SP,
                 para->getParD(0)->size_Mat_SP, para->getParD(1)->size_Mat_SP, para->getParD(0)->evenOrOdd,
                 para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF,
                 para->getParD(0)->K_CF, para->getParD(0)->omega, para->getParD(1)->omega,
                 para->getParD(0)->vis, para->getParD(0)->nx, para->getParD(0)->ny,
                 para->getParD(1)->nx, para->getParD(1)->ny, para->getParD(0)->numberofthreads,
                 para->getParD(0)->offCF);
             getLastCudaError("ScaleCF_RhoSq_comp_27 execution failed");


			//////////////////////////////////////////////////////////////////////////
           // ScaleCF_RhoSq_3rdMom_comp_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
										 //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
										 //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
										 //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
										 //para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
										 //para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
										 //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
										 //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
										 //para->getParD(0)->offCF);
           // getLastCudaError("ScaleCF_RhoSq_3rdMom_comp_27 execution failed");
            //////////////////////////////////////////////////////////////////////////
         //   ScaleCF_AA2016_comp_27( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
									//para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
									//para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
									//para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
									//para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
									//para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
									//para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
									//para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
									//para->getParD(0)->offCF);
         //   getLastCudaError("ScaleCF_AA2016_comp_27 execution failed");
            //////////////////////////////////////////////////////////////////////////


            //////////////////////////////////////////////////////////////////////////

            //ScaleCF_Fix_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
            //               para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offCF);
            //getLastCudaError("ScaleCF27 execution failed");

            //ScaleFC_Fix_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
            //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
            //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
            //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
            //               para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
            //               para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
            //               para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
            //               para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
            //               para->getParD(0)->offFC);
            //getLastCudaError("ScaleFC27 execution failed");
            //////////////////////////////////////////////////////////////////////////
          //  ScaleCF_NSPress_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
							   //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
							   //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
							   //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							   //para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
							   //para->getParD(0)->K_CF,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							   //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							   //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							   //para->getParD(0)->offCF);
          //  getLastCudaError("ScaleCF27 execution failed");

          //  ScaleFC_NSPress_27(para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0], 
							   //para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
							   //para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
							   //para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
							   //para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
							   //para->getParD(0)->K_FC,           para->getParD(0)->omega,           para->getParD(1)->omega, 
							   //para->getParD(0)->vis,            para->getParD(0)->nx,              para->getParD(0)->ny, 
							   //para->getParD(1)->nx,             para->getParD(1)->ny,              para->getParD(0)->numberofthreads,
							   //para->getParD(0)->offFC);
          //  getLastCudaError("ScaleFC27 execution failed");
            ////////////////////////////////////////////////////////////////////////
            if (para->getDiffOn()==true)
            {
               if (para->getDiffMod() == 7)
               {
                  //ScaleCFThS7(   para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
                  //               para->getParD(0)->d7.f[0],        para->getParD(1)->d7.f[0],                
                  //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
                  //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
                  //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                  //               para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
                  //               para->getParD(0)->K_CF,           
                  //               para->getParD(0)->vis,            para->getParD(1)->diffusivity,     para->getParD(0)->numberofthreads);
                  //getLastCudaError("ScaleCFTh7 execution failed");

                  //ScaleFCThS7(   para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],
                  //               para->getParD(0)->d7.f[0],        para->getParD(1)->d7.f[0],
                  //               para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
                  //               para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
                  //               para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                  //               para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
                  //               para->getParD(0)->K_FC,
                  //               para->getParD(0)->vis,            para->getParD(0)->diffusivity,     para->getParD(0)->numberofthreads);
                  //getLastCudaError("ScaleFCTh7 execution failed");
                  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  ScaleCFThSMG7( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
                                 para->getParD(0)->d7.f[0],        para->getParD(1)->d7.f[0],                
                                 para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
                                 para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
                                 para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                                 para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
                                 para->getParD(0)->K_CF,           
                                 para->getParD(0)->vis,            para->getParD(1)->diffusivity,     para->getParD(0)->numberofthreads,
                                 para->getParD(0)->offCF);
                  getLastCudaError("ScaleCFTh7 execution failed");

                  ScaleFCThSMG7( para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],
                                 para->getParD(0)->d7.f[0],        para->getParD(1)->d7.f[0],
                                 para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
                                 para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
                                 para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                                 para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
                                 para->getParD(0)->K_FC,
                                 para->getParD(0)->vis,            para->getParD(0)->diffusivity,     para->getParD(0)->numberofthreads,
                                 para->getParD(0)->offFC);
                  getLastCudaError("ScaleFCTh7 execution failed");
               } 
               else if (para->getDiffMod() == 27)
               {
                  ScaleCFThS27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],                
                                 para->getParD(0)->d27.f[0],       para->getParD(1)->d27.f[0],                
                                 para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP,
                                 para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP,
                                 para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                                 para->getParD(0)->intCF.ICellCFC, para->getParD(0)->intCF.ICellCFF, 
                                 para->getParD(0)->K_CF,           
								 para->getParD(0)->vis,            para->getParD(1)->diffusivity,     para->getParD(0)->numberofthreads,
								 para->getParD(0)->offCF);
                  getLastCudaError("ScaleCFTh27 execution failed");

                  ScaleFCThS27(  para->getParD(0)->d0SP.f[0],      para->getParD(1)->d0SP.f[0],
                                 para->getParD(0)->d27.f[0],       para->getParD(1)->d27.f[0],
                                 para->getParD(0)->neighborX_SP,   para->getParD(0)->neighborY_SP,    para->getParD(0)->neighborZ_SP, 
                                 para->getParD(1)->neighborX_SP,   para->getParD(1)->neighborY_SP,    para->getParD(1)->neighborZ_SP, 
                                 para->getParD(0)->size_Mat_SP,    para->getParD(1)->size_Mat_SP,     para->getParD(0)->evenOrOdd,
                                 para->getParD(0)->intFC.ICellFCC, para->getParD(0)->intFC.ICellFCF, 
                                 para->getParD(0)->K_FC,
								 para->getParD(0)->vis,            para->getParD(0)->diffusivity,     para->getParD(0)->numberofthreads,
								 para->getParD(0)->offFC);
                  getLastCudaError("ScaleFCTh27 execution failed");
               }
            } 
            //////////////////////////////////////////////////////////////////////////
         }
		 //output << "fertig \n";
      //}
      //+++++MPI++++++
      //if(numprocs>1){
      //   exchangeData();
      //}
	  //output << "vor CalcMed\n";
	  ////////////////////////////////////////////////////////////////////////////////




	  ////////////////////////////////////////////////////////////////////////////////
	  //Particles
	  ////////////////////////////////////////////////////////////////////////////////
	  if (para->getCalcParticle()) propagateParticles(para.get(), t);
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
				  para->cudaCopyFsForCheckPoint(lev);
			  }

			  output << "Dateien fuer CheckPoint schreiben t=" << t << "...";
			  restObj->safe(para.get());
			  rest->doCheckPoint(para->getFName(), t, para->getMyID());
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
										  para->getParD(lev)->RhoMP,		para->getParD(lev)->kMP,			para->getParD(lev)->numberOfPointskMP,
										  valuesPerClockCycle,				t_MP,								para->getParD(lev)->geoSP,
										  para->getParD(lev)->neighborX_SP, para->getParD(lev)->neighborY_SP,	para->getParD(lev)->neighborZ_SP,
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
				  para->cudaCopyMeasurePointsToHost(lev);
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
	  if (para->getDiffOn()==true) 
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
		  calcPlaneConc(para.get(), 0);
	  }
	  //////////////////////////////////////////////////////////////////////////////////




	  ////////////////////////////////////////////////////////////////////////////////
      // File IO
      ////////////////////////////////////////////////////////////////////////////////
      //comm.get()->startTimer();
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
		 //////////////////////////////////////////////////////////////////////////
		 //exchange data for valid post process
		 if (para->getNumprocs() > 1)
		 {
			 ////1D domain decomposition
			 //exchangePreCollDataGPU27(para, comm.get(), 0);
			 //3D domain decomposition
			 //output << "(print) start exchange Pre X (level 0) \n";
			 exchangePreCollDataXGPU27(para.get(), comm.get(), 0);
			 //output << "(print) end exchange Pre X (level 0) \n";
			 //output << "(print) start exchange Pre Y (level 0) \n";
			 exchangePreCollDataYGPU27(para.get(), comm.get(), 0);
			 //output << "(print) end exchange Pre Y (level 0) \n";
			 //output << "(print) start exchange Pre Z (level 0) \n";
			 exchangePreCollDataZGPU27(para.get(), comm.get(), 0);
			 //output << "(print) end exchange Pre Z (level 0) \n";
			 //////////////////////////////////////////////////////////////////////////
			 //3D domain decomposition convection diffusion
			 if (para->getDiffOn()==true)
			 {
				 exchangePreCollDataADXGPU27(para.get(), comm.get(), 0);
				 exchangePreCollDataADYGPU27(para.get(), comm.get(), 0);
				 exchangePreCollDataADZGPU27(para.get(), comm.get(), 0);
			 }
		 }
		 //////////////////////////////////////////////////////////////////////////

         if( para->getPrintFiles() )
         {
            output << "Dateien schreiben t=" << t << "...";
            ////////////////////////////////////////////////////////////////////////////////
            for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
            {
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

				   //überschreiben mit Wandknoten
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

			   para->cudaCopyPrint(lev);
			   if (para->getCalcMedian())
			   {
				   para->cudaCopyMedianPrint(lev);
			   }


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

				  para->cudaCopyConcDH(lev);
                  //cudaMemoryCopy(para->getParH(lev)->Conc, para->getParD(lev)->Conc,  para->getParH(lev)->mem_size_real_SP , cudaMemcpyDeviceToHost);
               }
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   ////print cp
			   //if ((para->getParH(lev)->cpTop.size() > 0) && (t > para->getTStartOut()))
			   //{
				  // printCpTopIntermediateStep(para, t, lev);
			   //}
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
			if (para->getCalc2ndOrderMoments())  calc2ndMoments(para.get());
			if (para->getCalc3rdOrderMoments())  calc3rdMoments(para.get());
			if (para->getCalcHighOrderMoments()) calcHigherOrderMoments(para.get());
			////////////////////////////////////////////////////////////////////////

			////////////////////////////////////////////////////////////////////////
			//calculate median on host
			////////////////////////////////////////////////////////////////////////
			if (para->getCalcMedian() && ((int)t > para->getTimeCalcMedStart()) && ((int)t <= para->getTimeCalcMedEnd()) && ((t%(unsigned int)para->getclockCycleForMP())==0))
			{
				unsigned int tdiff = t - t_prev;
				calcMedian(para.get(), tdiff);
			}
			////////////////////////////////////////////////////////////////////////
			writeTimestep(para.get(), t);
			////////////////////////////////////////////////////////////////////////
			//printDragLift(para, t);
			////////////////////////////////////////////////////////////////////////
			if (para->getCalcParticle()) copyAndPrintParticles(para.get(), t, false);
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

	////////////////////////////////////////////////////////////////////////////////
	//printDragLift(para);
	////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////
	if (para->getDiffOn()==true) printPlaneConc(para.get());
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



   //CudaFreeHostMemory
   checkCudaErrors( cudaEventDestroy(start_t));
   checkCudaErrors( cudaEventDestroy(stop_t));
   sdkDeleteTimer(&sdkTimer);
   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
   {
	  //para->cudaFreeFull(lev);
	  para->cudaFreeCoord(lev);
	  para->cudaFreeSP(lev);
	  if (para->getCalcMedian())
	  {
		  para->cudaFreeMedianSP(lev);
	  }
	  //para->cudaFreeVeloBC(lev);
	  //para->cudaFreeWallBC(lev);
	  //para->cudaFreeVeloBC(lev); 
	  para->cudaFreeInlet(lev);
	  para->cudaFreeOutlet(lev);
	  //para->cudaFreeGeomBC(lev);
	  //para->cudaFreePress(lev);
   }
   if (para->getMaxLevel()>1)
   {
      for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
      {
		 para->cudaFreeInterfaceCF(lev);
		 para->cudaFreeInterfaceFC(lev);
		 para->cudaFreeInterfaceOffCF(lev);
		 para->cudaFreeInterfaceOffFC(lev);
		 //para->cudaFreePressX1(lev);
      }
   }
   //para->cudaFreeVeloBC(0); //level = 0
   //para->cudaFreePressBC();
   //para->cudaFreeVeloPropeller(para->getFine());
   //para->cudaFreePressX0(para->getCoarse());

   //////////////////////////////////////////////////////////////////////////
   //Temp
   if (para->getDiffOn()==true)
   {
      for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
      {
         checkCudaErrors( cudaFreeHost(para->getParH(lev)->Conc_Full     ));
         checkCudaErrors( cudaFreeHost(para->getParH(lev)->Conc          ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->Temp.temp     ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->Temp.k        ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempVel.temp  ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempVel.velo  ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempVel.k     ));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempPress.temp));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempPress.velo));
		 checkCudaErrors( cudaFreeHost(para->getParH(lev)->TempPress.k   ));
      }
   }
   //////////////////////////////////////////////////////////////////////////


   //////////////////////////////////////////////////////////////////////////
   //free second order moments
   if (para->getCalc2ndOrderMoments())
   {
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   para->cudaFree2ndMoments(lev);
	   }
   }
   //////////////////////////////////////////////////////////////////////////
   //free third order moments
   if (para->getCalc3rdOrderMoments())
   {
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   para->cudaFree3rdMoments(lev);
	   }
   }
   //////////////////////////////////////////////////////////////////////////
   //free higher order moments
   if (para->getCalcHighOrderMoments())
   {
	   for (int lev=para->getCoarse(); lev <= para->getFine(); lev++)
	   {
		   para->cudaFreeHigherMoments(lev);
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
	   for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
	   {
		   //////////////////////////////////////////////////////////////////////////
		   for (unsigned int i=0; i < para->getNumberOfProcessNeighborsX(lev, "send"); i++)
		   {
			   para->cudaFreeProcessNeighborX(lev, i);
		   }
		   //////////////////////////////////////////////////////////////////////////
		   for (unsigned int i=0; i < para->getNumberOfProcessNeighborsY(lev, "send"); i++)
		   {
			   para->cudaFreeProcessNeighborY(lev, i);
		   }
		   //////////////////////////////////////////////////////////////////////////
		   for (unsigned int i=0; i < para->getNumberOfProcessNeighborsZ(lev, "send"); i++)
		   {
			   para->cudaFreeProcessNeighborZ(lev, i);
		   }
	   }
   }
   //////////////////////////////////////////////////////////////////////////
   //Normals
   if (para->getIsGeoNormal()){
	   for (int lev=para->getCoarse(); lev < para->getFine(); lev++)
	   {
		   para->cudaFreeGeomNormals(lev);
	   }
   }
   //////////////////////////////////////////////////////////////////////////
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
	pm[0] = new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level);
	pm[0]->setStartCoordinates(startX, startY, startZ);
	pm[0]->setEndCoordinates(endX, endY, endZ);
	pm[0]->setResistanceLBM();
	definePMarea(pm[0]);
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
	pm[1] = new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level);
	pm[1]->setStartCoordinates(startX, startY, startZ);
	pm[1]->setEndCoordinates(endX, endY, endZ);
	pm[1]->setResistanceLBM();
	definePMarea(pm[1]);
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
	pm[2] = new PorousMedia(porosity, geo, darcySI, forchheimerSI, dxLBM, dtLBM, level);
	pm[2]->setStartCoordinates(startX, startY, startZ);
	pm[2]->setEndCoordinates(endX, endY, endZ);
	pm[2]->setResistanceLBM();
	definePMarea(pm[2]);
	//////////////////////////////////////////////////////////////////////////

}

void Simulation::definePMarea(PorousMedia* pMedia)
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
	para->cudaCopySP(level);
	pMedia->setSizePM(counter);
	output << "definePMarea....cuda alloc PM \n";
	para->cudaAllocPorousMedia(pMedia, level);
	unsigned int *tpmArrayIDs = pMedia->getHostNodeIDsPM();
	
	output << "definePMarea....copy vector to array \n";
	for (unsigned int j = 0; j < pMedia->getSizePM(); j++)
	{
		tpmArrayIDs[j] = nodeIDsPorousMedia[j];
	}
	
	pMedia->setHostNodeIDsPM(tpmArrayIDs);
	output << "definePMarea....cuda copy PM \n";
	para->cudaCopyPorousMedia(pMedia, level);
}
