#include "Calculation/UpdateGrid27.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "Calculation/DragLift.h"
#include "Calculation/Cp.h"
//#include "Utilities/StringUtil.hpp"
//#include "Output/UnstructuredGridWriter.hpp"
#include "Communication/ExchangeData27.h"




void updateGrid27(Parameter* para, Communicator* comm, PorousMedia** pm, int level, int max_level, unsigned int t)
{
   if ( level == para->getFine() )
   {
      for (int internaltimestep=0;internaltimestep<2;internaltimestep++)
      {
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //if (t>para->getStartTurn()){
			 // para->setPhi((para->getPhi()+para->getParD(level)->deltaPhi));
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // KernelKum1hSP27(    para->getParD(level)->numberofthreads,       
				//				  para->getParD(level)->omega,
				//				  para->getParD(level)->deltaPhi,
				//				  para->getAngularVelocity(),
				//				  para->getParD(level)->geoSP, 
				//				  para->getParD(level)->neighborX_SP, 
				//				  para->getParD(level)->neighborY_SP, 
				//				  para->getParD(level)->neighborZ_SP,
				//				  para->getParD(level)->coordX_SP, 
				//				  para->getParD(level)->coordY_SP, 
				//				  para->getParD(level)->coordZ_SP,
				//				  para->getParD(level)->d0SP.f[0],    
				//				  para->getParD(level)->size_Mat_SP,  
				//				  para->getParD(level)->evenOrOdd); 
			 // getLastCudaError("KernelCasSPKum27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
				//				para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
				//				para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
				//				para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
				//				para->getPhi(),                        para->getAngularVelocity(),
				//				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
				//				para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
				//				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //    getLastCudaError("QVelDev27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //}
		  //else{
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // KernelKum1hSP27(    para->getParD(level)->numberofthreads,       
				//				  para->getParD(level)->omega,
				//				  (real)0.0,
				//				  (real)0.0,
				//				  para->getParD(level)->geoSP, 
				//				  para->getParD(level)->neighborX_SP, 
				//				  para->getParD(level)->neighborY_SP, 
				//				  para->getParD(level)->neighborZ_SP,
				//				  para->getParD(level)->coordX_SP, 
				//				  para->getParD(level)->coordY_SP, 
				//				  para->getParD(level)->coordZ_SP,
				//				  para->getParD(level)->d0SP.f[0],    
				//				  para->getParD(level)->size_Mat_SP,  
				//				  para->getParD(level)->evenOrOdd); 
			 // getLastCudaError("KernelCasSPKum27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
				//				para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
				//				para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
				//				para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
				//				para->getPhi(),                        (real)0.0,
				//				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
				//				para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
				//				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //    getLastCudaError("QVelDev27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //}
		  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelKumAA2016CompBulkSP27(para->getParD(level)->numberofthreads,       
				//						 para->getParD(level)->omega, 
				//						 para->getParD(level)->geoSP, 
				//						 para->getParD(level)->neighborX_SP, 
				//						 para->getParD(level)->neighborY_SP, 
				//						 para->getParD(level)->neighborZ_SP,
				//						 para->getParD(level)->d0SP.f[0],    
				//						 para->getParD(level)->size_Mat_SP,
				//						 para->getParD(level)->size_Array_SP,
				//						 level,
				//						 para->getForcesDev(),
				//						 para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelKumAA2016CompBulkSP27 execution failed");
			 //KernelKumAA2016CompSP27(para->getParD(level)->numberofthreads,       
				//				     para->getParD(level)->omega, 
				//				     para->getParD(level)->geoSP, 
				//				     para->getParD(level)->neighborX_SP, 
				//				     para->getParD(level)->neighborY_SP, 
				//				     para->getParD(level)->neighborZ_SP,
				//				     para->getParD(level)->d0SP.f[0],    
				//				     para->getParD(level)->size_Mat_SP,
				//				     level,
				//				     para->getForcesDev(),
				//				     para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelKumAA2016CompSP27 execution failed");
		 //KernelBGKPlusCompSP27(para->getParD(level)->numberofthreads,       
		 //					   para->getParD(level)->omega, 
		 //					   para->getParD(level)->geoSP, 
		 //					   para->getParD(level)->neighborX_SP, 
		 //					   para->getParD(level)->neighborY_SP, 
		 //					   para->getParD(level)->neighborZ_SP,
		 //					   para->getParD(level)->d0SP.f[0],    
		 //					   para->getParD(level)->size_Mat_SP,  
		 //					   para->getParD(level)->evenOrOdd); 
		 //getLastCudaError("KernelBGKPlusSP27 execution failed");
		 //printf("Level: %d \n", level);


		//KernelKumNewCompSP27(para->getParD(level)->numberofthreads,       
		//					  para->getParD(level)->omega, 
		//					  para->getParD(level)->geoSP, 
		//					  para->getParD(level)->neighborX_SP, 
		//					  para->getParD(level)->neighborY_SP, 
		//					  para->getParD(level)->neighborZ_SP,
		//					  para->getParD(level)->d0SP.f[0],    
		//					  para->getParD(level)->size_Mat_SP,  
		//					  para->getParD(level)->size_Array_SP,
		//					  level,
		//					  para->getForcesDev(),
		//					  para->getParD(level)->evenOrOdd); 
		//getLastCudaError("KernelCasSPKum27 execution failed");
 			//F3
			//KernelCumulantD3Q27F3(para->getParD(level)->numberofthreads,
			//					  para->getParD(level)->omega, 
			//					  para->getParD(level)->geoSP, 
			//					  para->getParD(level)->neighborX_SP, 
			//					  para->getParD(level)->neighborY_SP, 
			//					  para->getParD(level)->neighborZ_SP,
			//					  para->getParD(level)->d0SP.f[0],    
			//					  para->getParD(level)->g6.g[0],    
			//					  para->getParD(level)->size_Mat_SP,
			//					  level,
			//					  para->getForcesDev(),
			//					  para->getParD(level)->evenOrOdd);
			//getLastCudaError("KernelCumulantD3Q27F3 execution failed");

		  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //KernelKumCompSP27(  para->getParD(level)->numberofthreads,       
				//			  para->getParD(level)->omega, 
				//			  para->getParD(level)->geoSP, 
				//			  para->getParD(level)->neighborX_SP, 
				//			  para->getParD(level)->neighborY_SP, 
				//			  para->getParD(level)->neighborZ_SP,
				//			  para->getParD(level)->d0SP.f[0],    
				//			  para->getParD(level)->size_Mat_SP,  
				//			  para->getParD(level)->evenOrOdd); 
		  //getLastCudaError("KernelCasSPKum27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 		  //KernelBGKPlusSP27(para->getParD(level)->numberofthreads,       
 		 	//			    para->getParD(level)->omega, 
 		 	//			    para->getParD(level)->geoSP, 
 		 	//			    para->getParD(level)->neighborX_SP, 
 		 	//			    para->getParD(level)->neighborY_SP, 
 		 	//			    para->getParD(level)->neighborZ_SP,
 		 	//			    para->getParD(level)->d0SP.f[0],    
 		 	//			    para->getParD(level)->size_Mat_SP,  
 		 	//			    para->getParD(level)->evenOrOdd); 
 		  //getLastCudaError("KernelBGKPlusSP27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  		 // KernelMRTSP27(para->getParD(level)->numberofthreads,       
  			//			para->getParD(level)->omega, 
  			//			para->getParD(level)->geoSP, 
  			//			para->getParD(level)->neighborX_SP, 
  			//			para->getParD(level)->neighborY_SP, 
  			//			para->getParD(level)->neighborZ_SP,
  			//			para->getParD(level)->d0SP.f[0],    
  			//			para->getParD(level)->size_Mat_SP,  
  			//			para->getParD(level)->evenOrOdd); 
  		 //getLastCudaError("KernelMRT27 execution failed");
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 		 //KernelCascadeSP27(para->getParD(level)->numberofthreads,       
 			//			   para->getParD(level)->omega, 
 			//			   para->getParD(level)->geoSP, 
 			//			   para->getParD(level)->neighborX_SP, 
 			//			   para->getParD(level)->neighborY_SP, 
 			//			   para->getParD(level)->neighborZ_SP,
 			//			   para->getParD(level)->d0SP.f[0],    
 			//			   para->getParD(level)->size_Mat_SP,  
 			//			   para->getParD(level)->evenOrOdd); 
 		 // getLastCudaError("KernelCas27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //KernelKumNewSP27(   para->getParD(level)->numberofthreads,       
				//			  para->getParD(level)->omega, 
				//			  para->getParD(level)->geoSP, 
				//			  para->getParD(level)->neighborX_SP, 
				//			  para->getParD(level)->neighborY_SP, 
				//			  para->getParD(level)->neighborZ_SP,
				//			  para->getParD(level)->d0SP.f[0],    
				//			  para->getParD(level)->size_Mat_SP,  
				//			  para->getParD(level)->evenOrOdd); 
		  //getLastCudaError("KernelCasSPKum27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //KernelCasSPMSOHM27(  para->getParD(level)->numberofthreads,       
         //                     para->getParD(level)->omega, 
         //                     para->getParD(level)->geoSP, 
         //                     para->getParD(level)->neighborX_SP, 
         //                     para->getParD(level)->neighborY_SP, 
         //                     para->getParD(level)->neighborZ_SP,
         //                     para->getParD(level)->d0SP.f[0],    
         //                     para->getParD(level)->size_Mat_SP,  
         //                     para->getParD(level)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //KernelCasSPMS27(  para->getParD(level)->numberofthreads,       
         //                  para->getParD(level)->omega, 
         //                  para->getParD(level)->geoSP, 
         //                  para->getParD(level)->neighborX_SP, 
         //                  para->getParD(level)->neighborY_SP, 
         //                  para->getParD(level)->neighborZ_SP,
         //                  para->getParD(level)->d0SP.f[0],    
         //                  para->getParD(level)->size_Mat_SP,  
         //                  para->getParD(level)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //KernelCasSP27( para->getParD(level)->numberofthreads,       
         //               para->getParD(level)->omega, 
         //               para->getParD(level)->geoSP, 
         //               para->getParD(level)->neighborX_SP, 
         //               para->getParD(level)->neighborY_SP, 
         //               para->getParD(level)->neighborZ_SP,
         //               para->getParD(level)->d0SP.f[0],    
         //               para->getParD(level)->size_Mat_SP,  
         //               para->getParD(level)->evenOrOdd); 
         //getLastCudaError("KernelCasSP27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //porous media
		 if (para->getSimulatePorousMedia())
		 {
			 //////////////////////////////////////////////////////////////////////////
			 //Porous Media 0
			 KernelPMCumOneCompSP27(para->getParD(level)->numberofthreads,
									para->getParD(level)->omega,
									para->getParD(level)->neighborX_SP,
									para->getParD(level)->neighborY_SP,
									para->getParD(level)->neighborZ_SP,
									para->getParD(level)->d0SP.f[0],
									para->getParD(level)->size_Mat_SP,
									level,
									para->getForcesDev(),
									pm[0]->getPorosity(),
									pm[0]->getDarcyLBM(),
									pm[0]->getForchheimerLBM(),
									pm[0]->getSizePM(),
									pm[0]->getHostNodeIDsPM(),
									para->getParD(level)->evenOrOdd);
			 getLastCudaError("KernelPMCumOneCompSP27 execution failed");
			 //////////////////////////////////////////////////////////////////////////
			 //Porous Media 1
			 KernelPMCumOneCompSP27(para->getParD(level)->numberofthreads,
									para->getParD(level)->omega,
									para->getParD(level)->neighborX_SP,
									para->getParD(level)->neighborY_SP,
									para->getParD(level)->neighborZ_SP,
									para->getParD(level)->d0SP.f[0],
									para->getParD(level)->size_Mat_SP,
									level,
									para->getForcesDev(),
									pm[1]->getPorosity(),
									pm[1]->getDarcyLBM(),
									pm[1]->getForchheimerLBM(),
									pm[1]->getSizePM(),
									pm[1]->getHostNodeIDsPM(),
									para->getParD(level)->evenOrOdd);
			 getLastCudaError("KernelPMCumOneCompSP27 execution failed");
			 //////////////////////////////////////////////////////////////////////////
			 //Porous Media 2
			 KernelPMCumOneCompSP27(para->getParD(level)->numberofthreads,
									para->getParD(level)->omega,
									para->getParD(level)->neighborX_SP,
									para->getParD(level)->neighborY_SP,
									para->getParD(level)->neighborZ_SP,
									para->getParD(level)->d0SP.f[0],
									para->getParD(level)->size_Mat_SP,
									level,
									para->getForcesDev(),
									pm[2]->getPorosity(),
									pm[2]->getDarcyLBM(),
									pm[2]->getForchheimerLBM(),
									pm[2]->getSizePM(),
									pm[2]->getHostNodeIDsPM(),
									para->getParD(level)->evenOrOdd);
			 getLastCudaError("KernelPMCumOneCompSP27 execution failed");
		 }
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (para->getDiffOn())
         {
            if (para->getDiffMod() == 7)
            {
				//output << " Diff Mod 7\n";
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // incomp
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelADincomp7(para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
							   para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
							   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,  
							   para->getParD(level)->evenOrOdd); 
			   getLastCudaError("KernelADincomp7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion boundary condition
               QNoSlipADincompDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
									para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,  
									para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
									para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
									para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
									para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
               getLastCudaError("QNoSlipADincompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
			   if (t<15580)//(t>500000 && t<515580)//(t>300000 && t<315580)
			   {
                 QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
								   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
								   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
								   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
			   else
			   {
                 QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
								   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
								   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
								   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
           //    QADVeloIncompDev7(  para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   //para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
								   //para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   //para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	    para->getParD(level)->TempVel.kTemp,  
								   //para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	    para->getParD(level)->neighborY_SP, 
								   //para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
           //    getLastCudaError("QADVeloIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
               QADPressIncompDev7(   para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
									 para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
									 para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
									 para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp, para->getParD(level)->TempPress.kTemp,  
									 para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
									 para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
               getLastCudaError("QADPressIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			   

			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// comp
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion kernel
      //         KernelThS7( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
      //                     para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                     para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,  
      //                     para->getParD(level)->evenOrOdd); 
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         QADDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
      //                  para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,  
      //                  para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
      //                  para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
      //                  para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                  para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADVelDev7( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
      //                     para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
      //                     para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
      //                     para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,     para->getParD(level)->TempVel.kTemp,  
      //                     para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //                     para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADVelDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADPressDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
      //                       para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
      //                       para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
      //                       para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,   para->getParD(level)->TempPress.kTemp,  
      //                       para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //                       para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            } 
            else if (para->getDiffMod() == 27)
            {
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // incomp
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
      //         KernelADincomp27( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
						//		   para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
						//		   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],     para->getParD(level)->size_Mat_SP,  
						//		   para->getParD(level)->evenOrOdd); 
			   //getLastCudaError("KernelADincomp27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         QNoSlipADincompDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
						//			  para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,  
						//			  para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
						//			  para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
						//			  para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
						//			  para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QNoSlipADincompDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
			   //if (t>500000 && t<515580)//(t>300000 && t<315580)
			   //{
				   //QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				 para->getParD(level)->ny,
							//		    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		 para->getParD(level)->TempVel.tempPulse, 
							//		    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	 para->getParD(level)->TempVel.k,
							//		    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
							//		    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	 para->getParD(level)->neighborY_SP, 
							//		    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	 para->getParD(level)->evenOrOdd);
				   //getLastCudaError("QADVeloIncompDev27 execution failed");
			   //}
			   //else
			   //{
				  // QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//			    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//			    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	para->getParD(level)->TempVel.k,
						//			    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,  para->getParD(level)->TempVel.kTemp,  
						//			    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//			    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVeloIncompDev27 execution failed");
			   //}
      //     //    QADVeloIncompDev27( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//		     //para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//		     //para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
						//		     //para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
						//		     //para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//		     //para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //     //    getLastCudaError("QADVeloIncompDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + pressure boundary condition
           //    QADPressIncompDev27(   para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
									  //para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
									  //para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
									  //para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
									  //para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
									  //para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
           //    getLastCudaError("QADPressIncompDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// comp
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelThS27(para->getParD(level)->numberofthreads,
                           para->getParD(level)->diffusivity, 
                           para->getParD(level)->geoSP, 
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->d0SP.f[0],    
                           para->getParD(level)->d27.f[0],    
                           para->getParD(level)->size_Mat_SP,  
                           para->getParD(level)->evenOrOdd); 
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion boundary condition
               QADBBDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
                          para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,  
                          para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
                          para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                          para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
               getLastCudaError("QADBBDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //Street Manhattan - never use again, please
               QADDirichletDev27( para->getParD(level)->numberofthreads,      para->getParD(level)->nx,					para->getParD(level)->ny,
								  para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],			para->getParD(level)->TempVel.tempPulse,  
								  para->getParD(level)->diffusivity,          para->getParD(level)->concIndex,			para->getParD(level)->QGeom.q27[0], 
								  para->getParD(level)->numberOfPointsConc,   para->getParD(level)->numberOfPointsConc, para->getParD(level)->omega,
								  para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP,		para->getParD(level)->neighborZ_SP,
								  para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
               getLastCudaError("QADDirichletDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
			   //if (t<1000)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
			   //{
				  // QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//	   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.tempPulse, 
						//	   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
						//	   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,  
						//	   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//	   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
			   //else
			   //{
				  // QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//	   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//	   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
						//	   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,  
						//	   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//	   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         //QADPressDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
      //         //               para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
      //         //               para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
      //         //               para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
      //         //               para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //         //               para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         //getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
         }
		//////////////////////////////////////////////////////////////////////////////
		if (para->getNumprocs() > 1)
		{
			////1D domain decomposition
			//exchangePostCollDataGPU27(para, comm, level);
			//3D domain decomposition
			//printf("start exchange Post X (level: %d, myID: %d) \n", level, para->getMyID());
			exchangePostCollDataXGPU27(para, comm, level);
			//printf("end exchange Post X (level: %d, myID: %d) \n", level, para->getMyID());
			//printf("start exchange Post Y (level: %d, myID: %d) \n", level, para->getMyID());
			exchangePostCollDataYGPU27(para, comm, level);
			//printf("end exchange Post Y (level: %d, myID: %d) \n", level, para->getMyID());
			//printf("start exchange Post Z (level: %d, myID: %d) \n", level, para->getMyID());
			exchangePostCollDataZGPU27(para, comm, level);
			//printf("end exchange Post Z (level: %d, myID: %d) \n", level, para->getMyID());
			//////////////////////////////////////////////////////////////////////////
			//3D domain decomposition convection diffusion
			if (para->getDiffOn()==true)
			{
				exchangePostCollDataADXGPU27(para, comm, level);
				exchangePostCollDataADYGPU27(para, comm, level);
				exchangePostCollDataADZGPU27(para, comm, level);
			}
		}
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //QPressDevFixBackflow27( para->getParD(level)->numberofthreads,       RhoBCOutflowD,
         //                        para->getParD(level)->d0SP.f[0],    QoutflowD.k, kOutflowQ,             para->getParD(level)->omega,
         //                        para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
         //                        para->getParD(level)->size_Mat_SP,  para->getParD(level)->evenOrOdd);
         //getLastCudaError("QPressDev27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//QPressDev27_IntBB(  para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
		//					para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,       para->getParD(level)->QPress.q27[0], 
		//					para->getParD(level)->QPress.kQ,       para->getParD(level)->QPress.kQ,      para->getParD(level)->omega,
		//					para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
		//					para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QPressDev27_IntBB fine execution failed");
		 //QPressDevOld27( para->getParD(level)->numberofthreads, para->getParD(level)->QpressX1.RhoBC, 
			//			 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QpressX1.k,  
			//			 para->getParD(level)->QpressX1.kN,     para->getParD(level)->QpressX1.kQ,  para->getParD(level)->omega,
			//			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QPressDev27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//if (para->getParD(level)->kSlipQ > 0)
			//{
			//	QSlipDevComp27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
			//					para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
			//					para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//					para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			//	getLastCudaError("QSlipDev27 execution failed");
			//	//QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
			//	//			 para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
			//	//			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	//			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			//	//getLastCudaError("Slip27 execution failed");
			//	//QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QGeom.k,
			//	//			 para->getParD(level)->QGeom.q27[0],    para->getParD(level)->QGeom.kQ,       para->getParD(level)->omega,
			//	//			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	//			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			//	//getLastCudaError("Slip27 execution failed");
			//}
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 // printf("fein vor WallBC\n");
		 //QDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			//	 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
			//	 para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
			//	 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QDev27 execution failed");
		  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  // printf("fein vor WallBC\n");
       //   BBDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
       //            para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
				   //para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
       //            para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
       //            para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
       //   getLastCudaError("BBDev27 (Wall) execution failed");
		 //QDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			//	 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
			//	 para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
			//	 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QDev27 (Wall) execution failed");
			//QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			//	       para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
			//	       para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
			//	       para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	       para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			//getLastCudaError("QDevComp27 (Wall) execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		      //QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
								//para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
								//para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
								//para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
								//para->getPhi(),                        para->getAngularVelocity(),
								//para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
								//para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
								//para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		      //getLastCudaError("QVelDev27 execution failed");
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //???????????????????????????????
      //      QVelDev27(  para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
      //                  para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
      //                  para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
						//para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
      //                  para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                  para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
      //      getLastCudaError("QVelDev27 execution failed");
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//QVelDevComp27(  para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
			//				para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
			//				para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
			//				para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
			//				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			//getLastCudaError("QVelDevComp27 execution failed");

        //if (para->getParD(level)->kInflowQ > 0)
        //{
        //    QVelDevCompZeroPress27(para->getParD(level)->numberofthreads, para->getParD(level)->nx, para->getParD(level)->ny,
        //        para->getParD(level)->Qinflow.Vx, para->getParD(level)->Qinflow.Vy, para->getParD(level)->Qinflow.Vz,
        //        para->getParD(level)->d0SP.f[0], para->getParD(level)->Qinflow.k, para->getParD(level)->Qinflow.q27[0],
        //        para->getParD(level)->kInflowQ, para->getParD(level)->Qinflow.kArray, para->getParD(level)->omega,
        //        para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //        para->getParD(level)->size_Mat_SP, para->getParD(level)->evenOrOdd);
        //    getLastCudaError("QVelDevCompZeroPress27 execution failed");
        //}
			//QVelDevCompZeroPress27(  para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
			//						 para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
			//						 para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
			//						 para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
			//						 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//						 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			//getLastCudaError("QVelDevCompZeroPress27 execution failed");

		//  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////Drag and Lift Part I
		//DragLiftPostD27( para->getParD(level)->d0SP.f[0], 
		//				 para->getParD(level)->QGeom.k, 
		//				 para->getParD(level)->QGeom.q27[0],
		//				 para->getParD(level)->QGeom.kQ, 
		//				 para->getParD(level)->DragPostX,
		//				 para->getParD(level)->DragPostY,
		//				 para->getParD(level)->DragPostZ,
		//				 para->getParD(level)->neighborX_SP,
		//				 para->getParD(level)->neighborY_SP,
		//				 para->getParD(level)->neighborZ_SP,
		//				 para->getParD(level)->size_Mat_SP, 
		//				 para->getParD(level)->evenOrOdd,
		//				 para->getParD(level)->numberofthreads);
		//getLastCudaError("DragLift27 execution failed"); 
		//  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       //   BBDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
       //            para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
				   //para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
       //            para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
       //            para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
       //   getLastCudaError("BBDev27 (Wall) execution failed");
		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //Sphere
		  //QDev27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		  //		 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		  //		 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		  //		 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		  //		 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QDev27 (Geom) execution failed");

		 //if (para->getParD(level)->QGeom.kQ > 0)
		 //{
			//  QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		 // 				 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		 // 				 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		 // 				 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		 // 				 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			//  getLastCudaError("QDevComp27 (Geom) execution failed");
			//  //QDev3rdMomentsComp27(  para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		 // 	//						 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		 // 	//						 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		 // 	//						 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		 // 	//						 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			//  //getLastCudaError("QDev3rdMomentsComp27 (Geom) execution failed");
		 //  //   QVelDevCompZeroPress27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
			//		//				 para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
			//		//				 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
			//		//				 para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
			//		//				 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//		//				 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		 //  //   getLastCudaError("QVelDevCompZeroPress27 execution failed");
		 //}
		  //QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		  //			 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		  //			 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		  //			 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		  //			 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QDevComp27 (Geom) execution failed");

		  //QSlipGeomDevComp27(para->getParD(level)->numberofthreads,     para->getParD(level)->d0SP.f[0],           para->getParD(level)->QGeom.k,
				//			 para->getParD(level)->QGeom.q27[0],        para->getParD(level)->QGeom.kQ,            para->getParD(level)->omega,
				//			 para->getParD(level)->QGeomNormalX.q27[0], para->getParD(level)->QGeomNormalY.q27[0], para->getParD(level)->QGeomNormalZ.q27[0],
				//			 para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,        para->getParD(level)->neighborZ_SP,
				//			 para->getParD(level)->size_Mat_SP,         para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QSlipGeomDev27 execution failed");

		  //QSlipNormDevComp27(para->getParD(level)->numberofthreads,     para->getParD(level)->d0SP.f[0],           para->getParD(level)->QGeom.k,
				//			 para->getParD(level)->QGeom.q27[0],        para->getParD(level)->QGeom.kQ,            para->getParD(level)->omega,
				//			 para->getParD(level)->QGeomNormalX.q27[0], para->getParD(level)->QGeomNormalY.q27[0], para->getParD(level)->QGeomNormalZ.q27[0],
				//			 para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,        para->getParD(level)->neighborZ_SP,
				//			 para->getParD(level)->size_Mat_SP,         para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QSlipGeomDev27 execution failed");

		  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //Car
		  //QVelDev27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
		  //          para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
		  //          para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
				//	para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		  //          para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		  //          para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QVelDev27 execution failed");

		     // QVelDevComp27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
							//para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
							//para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
							//para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
							//para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
							//para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		     // getLastCudaError("QVelDevComp27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (para->getParD(level)->evenOrOdd==true)  para->getParD(level)->evenOrOdd=false;
         else                                        para->getParD(level)->evenOrOdd=true;
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 ////Drag and Lift Part II
		 //DragLiftPreD27(  para->getParD(level)->d0SP.f[0], 
			//			  para->getParD(level)->QGeom.k, 
			//			  para->getParD(level)->QGeom.q27[0],
			//			  para->getParD(level)->QGeom.kQ, 
			//			  para->getParD(level)->DragPreX,
			//			  para->getParD(level)->DragPreY,
			//			  para->getParD(level)->DragPreZ,
			//			  para->getParD(level)->neighborX_SP,
			//			  para->getParD(level)->neighborY_SP,
			//			  para->getParD(level)->neighborZ_SP,
			//			  para->getParD(level)->size_Mat_SP, 
			//			  para->getParD(level)->evenOrOdd,
			//			  para->getParD(level)->numberofthreads);
		 //getLastCudaError("DragLift27 execution failed"); 
		 //////////////////////////////////////////////////////////////////////////////////
		 ////Calculation of Drag and Lift
		 //////////////////////////////////////////////////////////////////////////////////
		 ////if(t>para->getStartTurn())calcDragLift(para, level);
		 //calcDragLift(para, level);
		 //////////////////////////////////////////////////////////////////////////////////

		 ////////////////////////////////////////////////////////////////////////////////
		  //QPressDevNEQ27( para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
		  //				  para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
		  //				  para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,    para->getParD(level)->omega,
		  //				  para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		  //				  para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QPressDev27 execution failed");
		  ////////////////////////////////////////////////////////////////////////////////

		 //////////////////////////////////////////////////////////////////////////////////
		 ////test for drag crisis
		 //if (para->getParD(level)->QGeom.kQ > 0)
		 //{
			// QDevCompHighNu27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			// 				  para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
			// 				  para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
			// 				  para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			// 				  para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			// getLastCudaError("QDevComp27 (Geom) execution failed");
		 //}
		 //////////////////////////////////////////////////////////////////////////////////

		 //////////////////////////////////////////////////////////////////////////////////
		 ////Calculation of cp
		 //////////////////////////////////////////////////////////////////////////////////
		 //if(t > para->getTStartOut())
		 //{
			// ////////////////////////////////////////////////////////////////////////////////
			// CalcCPtop27(para->getParD(level)->d0SP.f[0], 
			//			 para->getParD(level)->cpTopIndex, 
			//			 para->getParD(level)->numberOfPointsCpTop, 
			//			 para->getParD(level)->cpPressTop,
			//			 para->getParD(level)->neighborX_SP,
			//			 para->getParD(level)->neighborY_SP,
			//			 para->getParD(level)->neighborZ_SP,
			//			 para->getParD(level)->size_Mat_SP, 
			//			 para->getParD(level)->evenOrOdd,
			//			 para->getParD(level)->numberofthreads);
			// //////////////////////////////////////////////////////////////////////////////////
			// //CalcCPbottom27(para->getParD(level)->d0SP.f[0],
			//	//			para->getParD(level)->cpBottomIndex, 
			//	//			para->getParD(level)->numberOfPointsCpBottom, 
			//	//			para->getParD(level)->cpPressBottom,
			//	//			para->getParD(level)->neighborX_SP,
			//	//			para->getParD(level)->neighborY_SP,
			//	//			para->getParD(level)->neighborZ_SP,
			//	//			para->getParD(level)->size_Mat_SP, 
			//	//			para->getParD(level)->evenOrOdd,
			//	//			para->getParD(level)->numberofthreads);
			// //////////////////////////////////////////////////////////////////////////////////
			// //CalcCPbottom27(para->getParD(level)->d0SP.f[0],
			//	//			para->getParD(level)->cpBottom2Index, 
			//	//			para->getParD(level)->numberOfPointsCpBottom2, 
			//	//			para->getParD(level)->cpPressBottom2,
			//	//			para->getParD(level)->neighborX_SP,
			//	//			para->getParD(level)->neighborY_SP,
			//	//			para->getParD(level)->neighborZ_SP,
			//	//			para->getParD(level)->size_Mat_SP, 
			//	//			para->getParD(level)->evenOrOdd,
			//	//			para->getParD(level)->numberofthreads);
			// //////////////////////////////////////////////////////////////////////////////////
			// calcCp(para, level);
			// ////////////////////////////////////////////////////////////////////////////////
		 //}
		 //////////////////////////////////////////////////////////////////////////////////





   //      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 ////printf("fein vor Propeller\n");
		 //PropVelo(	para->getParD(level)->numberofthreads,   para->getParD(level)->neighborX_SP,    
			//		para->getParD(level)->neighborY_SP,	     para->getParD(level)->neighborZ_SP,
			//		para->getParD(level)->QPropeller.RhoBC,  para->getParD(level)->QPropeller.Vx,   
			//		para->getParD(level)->QPropeller.Vy,     para->getParD(level)->QPropeller.Vz,
			//		para->getParD(level)->QPropeller.k,      para->getParD(level)->QPropeller.kQ,     
			//		para->getParD(level)->size_Mat_SP,       para->getParD(level)->geoSP,
			//		para->getParD(level)->d0SP.f[0],         para->getParD(level)->evenOrOdd);
		 //getLastCudaError("PropVelo execution failed");
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		 if (para->getNumprocs() > 1)
		 {
			 ////1D domain decomposition
			 //exchangePreCollDataGPU27(para, comm, level);
			 //3D domain decomposition
			 //printf("start exchange Pre X (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataXGPU27(para, comm, level);
			 //printf("end exchange Pre X (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Pre Y (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataYGPU27(para, comm, level);
			 //printf("end exchange Pre Y (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Pre Z (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataZGPU27(para, comm, level);
			 //printf("end exchange Pre Z (level: %d, myID: %d) \n", level, para->getMyID());
			 //////////////////////////////////////////////////////////////////////////
			 //3D domain decomposition convection diffusion
			 if (para->getDiffOn()==true)
			 {
				 exchangePreCollDataADXGPU27(para, comm, level);
				 exchangePreCollDataADYGPU27(para, comm, level);
				 exchangePreCollDataADZGPU27(para, comm, level);
			 }
		 }

		 //////////////////////////////////////////////////////////////////////////////////
		 //// File IO TEST
		 //////////////////////////////////////////////////////////////////////////////////
		 ////comm->startTimer();
		 //if(para->getTOut()>0 && t%para->getTOut()==0 && t>para->getStartTurn())
		 //{
			// //////////////////////////////////////////////////////////////////////////
			// //Timer SDK
			// //////////////////////////////////////////////////////////////////////////
			// if( para->getPrintFiles() )
			// {
			//	CalcMacSP27(para->getParD(level)->vx_SP,       
			//				 para->getParD(level)->vy_SP,        
			//				 para->getParD(level)->vz_SP,        
			//				 para->getParD(level)->rho_SP, 
			//				 para->getParD(level)->press_SP, 
			//				 para->getParD(level)->geoSP,       
			//				 para->getParD(level)->neighborX_SP, 
			//				 para->getParD(level)->neighborY_SP, 
			//				 para->getParD(level)->neighborZ_SP,
			//				 para->getParD(level)->size_Mat_SP, 
			//				 para->getParD(level)->numberofthreads,       
			//				 para->getParD(level)->d0SP.f[0],    
			//				 para->getParD(level)->evenOrOdd);
			//	getLastCudaError("CalcMacSP27 execution failed"); 

			//	para->cudaCopyPrint(level);

			//	std::string ffname_bin = para->getFName()+"_bin_"+StringUtil::toString<int>(level)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+"_"+StringUtil::toString<int>(internaltimestep)+".vtk";
			//	std::string ffname_bin_Points = para->getFName()+"_Points_"+"_bin_"+StringUtil::toString<int>(level)+"_"+StringUtil::toString<int>(para->getMyID())+"_"+StringUtil::toString<int>(t)+"_"+StringUtil::toString<int>(internaltimestep)+".vtk";
			//	UnstrucuredGridWriter::writeUnstrucuredGridEff(para, level, ffname_bin, ffname_bin_Points);
			//	////////////////////////////////////////////////////////////////////////////////
			//}
		 //}
      }
   }
   else
   {
      for (int internaltimestep=0;internaltimestep<2;internaltimestep++)
      {
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         updateGrid27(para, comm, pm, level+1, max_level, t);
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //if (t>para->getStartTurn()){
			 // //para->setPhi((para->getPhi()+para->getParD(level)->deltaPhi));
			 // KernelKum1hSP27(    para->getParD(level)->numberofthreads,       
				//				  para->getParD(level)->omega,
				//				  para->getParD(level)->deltaPhi,
				//				  para->getAngularVelocity(),
				//				  para->getParD(level)->geoSP, 
				//				  para->getParD(level)->neighborX_SP, 
				//				  para->getParD(level)->neighborY_SP, 
				//				  para->getParD(level)->neighborZ_SP,
				//				  para->getParD(level)->coordX_SP, 
				//				  para->getParD(level)->coordY_SP, 
				//				  para->getParD(level)->coordZ_SP,
				//				  para->getParD(level)->d0SP.f[0],    
				//				  para->getParD(level)->size_Mat_SP,  
				//				  para->getParD(level)->evenOrOdd); 
			 // getLastCudaError("KernelCasSPKum27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
				//				para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
				//				para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
				//				para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
				//				para->getPhi(),                        para->getAngularVelocity(),
				//				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
				//				para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
				//				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //    getLastCudaError("QVelDev27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //}
		  //else{
			 // KernelKum1hSP27(    para->getParD(level)->numberofthreads,       
				//				  para->getParD(level)->omega,
				//				  (real)0.0,
				//				  (real)0.0,
				//				  para->getParD(level)->geoSP, 
				//				  para->getParD(level)->neighborX_SP, 
				//				  para->getParD(level)->neighborY_SP, 
				//				  para->getParD(level)->neighborZ_SP,
				//				  para->getParD(level)->coordX_SP, 
				//				  para->getParD(level)->coordY_SP, 
				//				  para->getParD(level)->coordZ_SP,
				//				  para->getParD(level)->d0SP.f[0],    
				//				  para->getParD(level)->size_Mat_SP,  
				//				  para->getParD(level)->evenOrOdd); 
			 // getLastCudaError("KernelCasSPKum27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
				//				para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
				//				para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
				//				para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
				//				para->getPhi(),                        (real)0.0,
				//				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
				//				para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
				//				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		  //    getLastCudaError("QVelDev27 execution failed");
			 // //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		  //}
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //Sponge Test
		 //if (level == 7)
		 //{
			// KernelKumNewCompSpongeSP27(para->getParD(level)->numberofthreads,       
			//							para->getParD(level)->omega, 
			//							para->getParD(level)->geoSP, 
			//							para->getParD(level)->neighborX_SP, 
			//							para->getParD(level)->neighborY_SP, 
			//							para->getParD(level)->neighborZ_SP,
			//							para->getParD(level)->coordX_SP,
			//							para->getParD(level)->coordY_SP,
			//							para->getParD(level)->coordZ_SP,
			//							para->getParD(level)->d0SP.f[0],    
			//							para->getParD(level)->size_Mat_SP,  
			//							para->getParD(level)->evenOrOdd); 
			// getLastCudaError("KernelCasSPKum27 execution failed");
		 //}
		 //else
		 //{
			// KernelKumNewCompSP27(para->getParD(level)->numberofthreads,       
			//					  para->getParD(level)->omega, 
			//					  para->getParD(level)->geoSP, 
			//					  para->getParD(level)->neighborX_SP, 
			//					  para->getParD(level)->neighborY_SP, 
			//					  para->getParD(level)->neighborZ_SP,
			//					  para->getParD(level)->d0SP.f[0],    
			//					  para->getParD(level)->size_Mat_SP,  
			//					  para->getParD(level)->evenOrOdd); 
			// getLastCudaError("KernelCasSPKum27 execution failed");
		 //}
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelKumAA2016CompBulkSP27(para->getParD(level)->numberofthreads,       
				//						 para->getParD(level)->omega, 
				//						 para->getParD(level)->geoSP, 
				//						 para->getParD(level)->neighborX_SP, 
				//						 para->getParD(level)->neighborY_SP, 
				//						 para->getParD(level)->neighborZ_SP,
				//						 para->getParD(level)->d0SP.f[0],    
				//						 para->getParD(level)->size_Mat_SP,
				//						 para->getParD(level)->size_Array_SP,
				//						 level,
				//						 para->getForcesDev(),
				//						 para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelKumAA2016CompBulkSP27 execution failed");
			 //KernelKumAA2016CompSP27(para->getParD(level)->numberofthreads,       
				//				     para->getParD(level)->omega, 
				//				     para->getParD(level)->geoSP, 
				//				     para->getParD(level)->neighborX_SP, 
				//				     para->getParD(level)->neighborY_SP, 
				//				     para->getParD(level)->neighborZ_SP,
				//				     para->getParD(level)->d0SP.f[0],    
				//				     para->getParD(level)->size_Mat_SP,
				//				     level,
				//				     para->getForcesDev(),
				//				     para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelKumAA2016CompSP27 execution failed");
			 //KernelBGKPlusCompSP27(para->getParD(level)->numberofthreads,       
			 //					   para->getParD(level)->omega, 
			 //					   para->getParD(level)->geoSP, 
			 //					   para->getParD(level)->neighborX_SP, 
			 //					   para->getParD(level)->neighborY_SP, 
			 //					   para->getParD(level)->neighborZ_SP,
			 //					   para->getParD(level)->d0SP.f[0],    
			 //					   para->getParD(level)->size_Mat_SP,  
			 //					   para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelBGKPlusSP27 execution failed");
			 //printf("Level: %d \n", level);
			//KernelKumNewCompSP27(para->getParD(level)->numberofthreads,       
			//					  para->getParD(level)->omega, 
			//					  para->getParD(level)->geoSP, 
			//					  para->getParD(level)->neighborX_SP, 
			//					  para->getParD(level)->neighborY_SP, 
			//					  para->getParD(level)->neighborZ_SP,
			//					  para->getParD(level)->d0SP.f[0],    
			//					  para->getParD(level)->size_Mat_SP,
			//					  para->getParD(level)->size_Array_SP,
			//					  level,
			//					  para->getForcesDev(),
			//					  para->getParD(level)->evenOrOdd); 
			//getLastCudaError("KernelCasSPKum27 execution failed");
			//KernelCumulantD3Q27F3(para->getParD(level)->numberofthreads,
			//					  para->getParD(level)->omega, 
			//					  para->getParD(level)->geoSP, 
			//					  para->getParD(level)->neighborX_SP, 
			//					  para->getParD(level)->neighborY_SP, 
			//					  para->getParD(level)->neighborZ_SP,
			//					  para->getParD(level)->d0SP.f[0],    
			//					  para->getParD(level)->g6.g[0],    
			//					  para->getParD(level)->size_Mat_SP,
			//					  level,
			//					  para->getForcesDev(),
			//					  para->getParD(level)->evenOrOdd);
			//getLastCudaError("KernelCumulantD3Q27F3 execution failed");
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelKumCompSP27(  para->getParD(level)->numberofthreads,       
			 //			  para->getParD(level)->omega, 
			 //			  para->getParD(level)->geoSP, 
			 //			  para->getParD(level)->neighborX_SP, 
			 //			  para->getParD(level)->neighborY_SP, 
			 //			  para->getParD(level)->neighborZ_SP,
			 //			  para->getParD(level)->d0SP.f[0],    
			 //			  para->getParD(level)->size_Mat_SP,  
			 //			  para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelCasSPKum27 execution failed");
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelBGKPlusSP27(para->getParD(level)->numberofthreads,       
			 //			    para->getParD(level)->omega, 
			 //			    para->getParD(level)->geoSP, 
			 //			    para->getParD(level)->neighborX_SP, 
			 //			    para->getParD(level)->neighborY_SP, 
			 //			    para->getParD(level)->neighborZ_SP,
			 //			    para->getParD(level)->d0SP.f[0],    
			 //			    para->getParD(level)->size_Mat_SP,  
			 //			    para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelBGKPlusSP27 execution failed");
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 // KernelMRTSP27(para->getParD(level)->numberofthreads,       
			 //			para->getParD(level)->omega, 
			 //			para->getParD(level)->geoSP, 
			 //			para->getParD(level)->neighborX_SP, 
			 //			para->getParD(level)->neighborY_SP, 
			 //			para->getParD(level)->neighborZ_SP,
			 //			para->getParD(level)->d0SP.f[0],    
			 //			para->getParD(level)->size_Mat_SP,  
			 //			para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelMRT27 execution failed");
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelCascadeSP27(para->getParD(level)->numberofthreads,       
			 //			   para->getParD(level)->omega, 
			 //			   para->getParD(level)->geoSP, 
			 //			   para->getParD(level)->neighborX_SP, 
			 //			   para->getParD(level)->neighborY_SP, 
			 //			   para->getParD(level)->neighborZ_SP,
			 //			   para->getParD(level)->d0SP.f[0],    
			 //			   para->getParD(level)->size_Mat_SP,  
			 //			   para->getParD(level)->evenOrOdd); 
			 // getLastCudaError("KernelCas27 execution failed");
			 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelKumNewSP27(    para->getParD(level)->numberofthreads,       
			 //				  para->getParD(level)->omega, 
			 //				  para->getParD(level)->geoSP, 
			 //				  para->getParD(level)->neighborX_SP, 
			 //				  para->getParD(level)->neighborY_SP, 
			 //				  para->getParD(level)->neighborZ_SP,
			 //				  para->getParD(level)->d0SP.f[0],    
			 //				  para->getParD(level)->size_Mat_SP,  
			 //				  para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelCasSPKum27 execution failed");
			 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelCasSPMSOHM27(  para->getParD(level)->numberofthreads,       
			 //                     para->getParD(level)->omega, 
			 //                     para->getParD(level)->geoSP, 
			 //                     para->getParD(level)->neighborX_SP, 
			 //                     para->getParD(level)->neighborY_SP, 
			 //                     para->getParD(level)->neighborZ_SP,
			 //                     para->getParD(level)->d0SP.f[0],    
			 //                     para->getParD(level)->size_Mat_SP,  
			 //                     para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelCasSP27 execution failed");
			 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelCasSPMS27(     para->getParD(level)->numberofthreads,       
			 //                     para->getParD(level)->omega, 
			 //                     para->getParD(level)->geoSP, 
			 //                     para->getParD(level)->neighborX_SP, 
			 //                     para->getParD(level)->neighborY_SP, 
			 //                     para->getParD(level)->neighborZ_SP,
			 //                     para->getParD(level)->d0SP.f[0],    
			 //                     para->getParD(level)->size_Mat_SP,  
			 //                     para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelCasSP27 execution failed");
			 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			 //KernelCasSP27( para->getParD(level)->numberofthreads,       
			 //               para->getParD(level)->omega, 
			 //               para->getParD(level)->geoSP, 
			 //               para->getParD(level)->neighborX_SP, 
			 //               para->getParD(level)->neighborY_SP, 
			 //               para->getParD(level)->neighborZ_SP,
			 //               para->getParD(level)->d0SP.f[0],    
			 //               para->getParD(level)->size_Mat_SP,  
			 //               para->getParD(level)->evenOrOdd); 
			 //getLastCudaError("KernelCasSP27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //if (para->getDiffOn()==true)
		 //{
			// if (para->getDiffMod() == 7)
			// {
			//	 KernelThS7( para->getParD(level)->numberofthreads,
			//		 para->getParD(level)->diffusivity, 
			//		 para->getParD(level)->geoSP, 
			//		 para->getParD(level)->neighborX_SP, 
			//		 para->getParD(level)->neighborY_SP, 
			//		 para->getParD(level)->neighborZ_SP,
			//		 para->getParD(level)->d0SP.f[0],    
			//		 para->getParD(level)->d7.f[0],    
			//		 para->getParD(level)->size_Mat_SP,  
			//		 para->getParD(level)->evenOrOdd); 
			//	 getLastCudaError("KernelThS7 execution failed");
			// } 
			// else if (para->getDiffMod() == 27)
			// {
			//	 KernelThS27(para->getParD(level)->numberofthreads,
			//		 para->getParD(level)->diffusivity, 
			//		 para->getParD(level)->geoSP, 
			//		 para->getParD(level)->neighborX_SP, 
			//		 para->getParD(level)->neighborY_SP, 
			//		 para->getParD(level)->neighborZ_SP,
			//		 para->getParD(level)->d0SP.f[0],    
			//		 para->getParD(level)->d27.f[0],    
			//		 para->getParD(level)->size_Mat_SP,  
			//		 para->getParD(level)->evenOrOdd); 
			//	 getLastCudaError("KernelThS27 execution failed");
			// }
		 //}
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (para->getDiffOn()==true)
         {
            if (para->getDiffMod() == 7)
            {
				//output << " Diff Mod 7\n";
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // incomp
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelADincomp7(para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
							   para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
							   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,  
							   para->getParD(level)->evenOrOdd); 
			   getLastCudaError("KernelADincomp7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion boundary condition
               QNoSlipADincompDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
									para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,  
									para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
									para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
									para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
									para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
               getLastCudaError("QNoSlipADincompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
			   if (t<15580)//(t>500000 && t<515580)//(t>300000 && t<315580)
			   {
                 QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
								   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
								   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
								   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
			   else
			   {
                 QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
								   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
								   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
								   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                 getLastCudaError("QADVeloIncompDev7 execution failed");

			   }
           //    QADVeloIncompDev7(  para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
								   //para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
								   //para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
								   //para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	    para->getParD(level)->TempVel.kTemp,  
								   //para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	    para->getParD(level)->neighborY_SP, 
								   //para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
           //    getLastCudaError("QADVeloIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + velocity boundary condition
               QADPressIncompDev7(   para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
									 para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
									 para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
									 para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp, para->getParD(level)->TempPress.kTemp,  
									 para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
									 para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
               getLastCudaError("QADPressIncompDev7 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			   

			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// comp
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion kernel
      //         KernelThS7( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
      //                     para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                     para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,  
      //                     para->getParD(level)->evenOrOdd); 
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         QADDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
      //                  para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,  
      //                  para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
      //                  para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
      //                  para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                  para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADVelDev7( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
      //                     para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
      //                     para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
      //                     para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,     para->getParD(level)->TempVel.kTemp,  
      //                     para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //                     para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADVelDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         QADPressDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
      //                       para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
      //                       para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
      //                       para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,   para->getParD(level)->TempPress.kTemp,  
      //                       para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //                       para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            } 
            else if (para->getDiffMod() == 27)
            {
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   // incomp
			   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
      //         KernelADincomp27( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP, 
						//		   para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
						//		   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],     para->getParD(level)->size_Mat_SP,  
						//		   para->getParD(level)->evenOrOdd); 
			   //getLastCudaError("KernelADincomp27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion boundary condition
      //         QNoSlipADincompDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
						//			  para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,  
						//			  para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
						//			  para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
						//			  para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
						//			  para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
      //         getLastCudaError("QNoSlipADincompDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
			   //if (t>500000 && t<515580)//(t>300000 && t<315580)
			   //{
				   //QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				 para->getParD(level)->ny,
							//		    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		 para->getParD(level)->TempVel.tempPulse, 
							//		    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	 para->getParD(level)->TempVel.k,
							//		    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
							//		    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	 para->getParD(level)->neighborY_SP, 
							//		    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	 para->getParD(level)->evenOrOdd);
				   //getLastCudaError("QADVeloIncompDev27 execution failed");
			   //}
			   //else
			   //{
				  // QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//			    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//			    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	para->getParD(level)->TempVel.k,
						//			    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,  para->getParD(level)->TempVel.kTemp,  
						//			    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//			    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVeloIncompDev27 execution failed");
			   //}
      //     //    QADVeloIncompDev27( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//		     //para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//		     //para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
						//		     //para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
						//		     //para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//		     //para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //     //    getLastCudaError("QADVeloIncompDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion + pressure boundary condition
           //    QADPressIncompDev27(   para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
									  //para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
									  //para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
									  //para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
									  //para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
									  //para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
           //    getLastCudaError("QADPressIncompDev27 execution failed");
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



			   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			   //// comp
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion kernel
               KernelThS27(para->getParD(level)->numberofthreads,
                           para->getParD(level)->diffusivity, 
                           para->getParD(level)->geoSP, 
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->d0SP.f[0],    
                           para->getParD(level)->d27.f[0],    
                           para->getParD(level)->size_Mat_SP,  
                           para->getParD(level)->evenOrOdd); 
               ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
               //advection diffusion boundary condition
               QADBBDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
                          para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,  
                          para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0], 
                          para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                          para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
               getLastCudaError("QADBBDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
			   //if (t<1000)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
			   //{
				  // QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//	   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.tempPulse, 
						//	   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
						//	   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,  
						//	   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//	   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
			   //else
			   //{
				  // QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
						//	   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
						//	   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
						//	   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,  
						//	   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
						//	   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
				  // getLastCudaError("QADVelDev27 execution failed");
			   //}
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //         //advection diffusion + velocity boundary condition
      //         //QADPressDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
      //         //               para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
      //         //               para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
      //         //               para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
      //         //               para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
      //         //               para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
      //         //getLastCudaError("QADPressDev27 execution failed");
      //         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            }
         }
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 if (para->getNumprocs() > 1)
		 {
			 ////1D domain decomposition
			 //exchangePostCollDataGPU27(para, comm, level);
			 //3D domain decomposition
			 //printf("start exchange Post X (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePostCollDataXGPU27(para, comm, level);
			 //printf("end exchange Post X (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Post Y (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePostCollDataYGPU27(para, comm, level);
			 //printf("end exchange Post Y (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Post Z (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePostCollDataZGPU27(para, comm, level);
			 //printf("end exchange Post Z (level: %d, myID: %d) \n", level, para->getMyID());
			 //////////////////////////////////////////////////////////////////////////
			 //3D domain decomposition convection diffusion
			 if (para->getDiffOn()==true)
			 {
				 exchangePostCollDataADXGPU27(para, comm, level);
				 exchangePostCollDataADYGPU27(para, comm, level);
				 exchangePostCollDataADZGPU27(para, comm, level);
			 }
		 }
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//QPressDev27_IntBB(  para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
		//					para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,       para->getParD(level)->QPress.q27[0], 
		//					para->getParD(level)->QPress.kQ,       para->getParD(level)->QPress.kQ,      para->getParD(level)->omega,
		//					para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
		//					para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QPressDev27_IntBB fine execution failed");
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //QPressDevOld27( para->getParD(level)->numberofthreads, para->getParD(level)->QpressX1.RhoBC, 
			//			 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QpressX1.k,  
			//			 para->getParD(level)->QpressX1.kN,     para->getParD(level)->QpressX1.kQ,  para->getParD(level)->omega,
			//			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QPressDev27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //if (para->getParD(level)->kSlipQ > 0)
		 //{
			// QSlipDevComp27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
			//	 para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
			//	 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			// getLastCudaError("QSlipDev27 execution failed");
			// //QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
			// //			 para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
			// //			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			// //			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			// //getLastCudaError("Slip27 execution failed");
			// //QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QGeom.k,
			// //			 para->getParD(level)->QGeom.q27[0],    para->getParD(level)->QGeom.kQ,       para->getParD(level)->omega,
			// //			 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			// //			 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			// //getLastCudaError("Slip27 execution failed");
		 //}
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //if (para->getParD(level)->kQ > 0)
		 //{
			// QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			//		    para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
			//		    para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
			//			para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//			para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			//getLastCudaError("QDevComp27 (Wall) execution failed");
		 //}
		 //QDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
			//	 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
			//	 para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
			//	 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//	 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QDev27 execution failed");
		 if (para->getParD(level)->kInflowQ > 0)
		 {
		      QVelDevCompZeroPress27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,             para->getParD(level)->ny,
									 para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,     para->getParD(level)->Qinflow.Vz,
									 para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,      para->getParD(level)->Qinflow.q27[0], 
									 para->getParD(level)->kInflowQ,        para->getParD(level)->Qinflow.kArray, para->getParD(level)->omega,
									 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
									 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		      getLastCudaError("QVelDevCompZeroPress27 execution failed");
		 }

		 //if (para->getParD(level)->QGeom.kQ > 0)
		 //{
			//  QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		 // 				 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		 // 				 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		 // 				 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		 // 				 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
			//  getLastCudaError("QDevComp27 (Geom) execution failed");
		 //     //QVelDevCompZeroPress27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
			//					//	 para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
			//					//	 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
			//					//	 para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
			//					//	 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
			//					//	 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		 //     //getLastCudaError("QVelDevCompZeroPress27 execution failed");
		 //}
		  //QDev27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		  //		 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		  //		 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		  //		 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		  //		 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		  //getLastCudaError("QDev27 (Geom) execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		      //QVelDevice1h27(   para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
								//para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
								//para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
								//para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
								//para->getPhi(),                        para->getAngularVelocity(),
								//para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
								//para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
								//para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		      //getLastCudaError("QVelDev27 execution failed");
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //??????????????????????????????????
      //      QVelDev27(  para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
      //                  para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
      //                  para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
						//para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
      //                  para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
      //                  para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
      //      getLastCudaError("QVelDev27 execution failed");

		     // QVelDevComp27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
							//para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
							//para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
							//para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
							//para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
							//para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		     // getLastCudaError("QVelDevComp27 execution failed");

         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         if (para->getParD(level)->evenOrOdd==true)  para->getParD(level)->evenOrOdd=false;
         else                                        para->getParD(level)->evenOrOdd=true;
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		 //////////////////////////////////////////////////////////////////////////////////
		 //press NEQ comp
		 //QPressDevNEQ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
		 //				para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
		 //				para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,    para->getParD(level)->omega,
		 //				para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		 //				para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		 //getLastCudaError("QPressDevNEQ27 execution failed");
		 ////////////////////////////////////////////////////////////////////////////////
		 //press EQ comp
		 //if (para->getParD(level)->QPress.kQ > 0)
		 //{
			// QPressDevEQZ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
		 //					para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
		 //					para->getParD(level)->QPress.kN,       para->getParD(level)->kDistTestRE.f[0],       
			//				para->getParD(level)->QPress.kQ,       para->getParD(level)->omega,
		 //					para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		 //					para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
			// getLastCudaError("QPressDevEQZ27 execution failed");
		 //}
		 ////////////////////////////////////////////////////////////////////////////////



		 //Spiegelbenchmark
		 //////////////////////////////////////////////////////////////////////////////////
		 ////Calculation of cp
		 //////////////////////////////////////////////////////////////////////////////////
		 //if (/*(para->getParD(level)->numberOfPointsCpTop > 0)*/ (level == 5) && (t > para->getTStartOut()))
		 //{
			// ////////////////////////////////////////////////////////////////////////////////
			// //Level 7
			// CalcCPtop27(para->getParD(7)->d0SP.f[0],
			//			 para->getParD(7)->cpTopIndex,
			//			 para->getParD(7)->numberOfPointsCpTop,
			//			 para->getParD(7)->cpPressTop,
			//			 para->getParD(7)->neighborX_SP,
			//			 para->getParD(7)->neighborY_SP,
			//			 para->getParD(7)->neighborZ_SP,
			//			 para->getParD(7)->size_Mat_SP,
			//			 para->getParD(7)->evenOrOdd,
			//			 para->getParD(7)->numberofthreads);
			// //////////////////////////////////////////////////////////////////////////////////
			// calcPressForMirror(para, 7);
			// ////////////////////////////////////////////////////////////////////////////////
			// //Level 8
			// CalcCPtop27(para->getParD(8)->d0SP.f[0],
			//			 para->getParD(8)->cpTopIndex,
			//			 para->getParD(8)->numberOfPointsCpTop,
			//			 para->getParD(8)->cpPressTop,
			//			 para->getParD(8)->neighborX_SP,
			//			 para->getParD(8)->neighborY_SP,
			//			 para->getParD(8)->neighborZ_SP,
			//			 para->getParD(8)->size_Mat_SP,
			//			 para->getParD(8)->evenOrOdd,
			//			 para->getParD(8)->numberofthreads);
			// //////////////////////////////////////////////////////////////////////////////////
			// calcPressForMirror(para, 8);
			// ////////////////////////////////////////////////////////////////////////////////
			// //print press mirror
			// printScalars(para, true); //true for binary
			// ////////////////////////////////////////////////////////////////////////////////
		 //}
		 //////////////////////////////////////////////////////////////////////////////////




		 ////////////////////////////////////////////////////////////////////////////////
		 //fine to coarse interpolation
			 //ScaleFC_comp_D3Q27F3(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],     para->getParD(level)->g6.g[0],
				//			      para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
				//			      para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
				//			      para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
				//			      para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
				//			      para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
				//			      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
				//			      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
				//			      para->getParD(level)->offFC);
    //        getLastCudaError("ScaleFC_comp_D3Q27F3 execution failed");
			////////////////////////////////////////////////////////////////////////////////
			//ScaleFC_0817_comp_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
							     //para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
							     //para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
							     //para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
							     //para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
							     //para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
							     //para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
							     //para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
							     //para->getParD(level)->offFC);
            //getLastCudaError("ScaleFC_0817_comp_27 execution failed");
         //ScaleFC_RhoSq_3rdMom_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
									//	para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
									//	para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
									//	para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
									//	para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
									//	para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
									//	para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
									//	para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
									//	para->getParD(level)->offFC);
         //getLastCudaError("ScaleFC_RhoSq_3rdMom_comp_27 execution failed");
		 ////////////////////////////////////////////////////////////////////////////////
         ScaleFC_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
								para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
								para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
								para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
								para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								para->getParD(level)->offFC);
         getLastCudaError("ScaleFC27_RhoSq_comp execution failed");
		 ////////////////////////////////////////////////////////////////////////////////
		 //data exchange
		 if (para->getNumprocs() > 1)
		 {
			 ////1D domain decomposition
			 //exchangePreCollDataGPU27(para, comm, level);
			 //3D domain decomposition
			 //printf("start exchange Pre X (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataXGPU27(para, comm, level);
			 //printf("end exchange Pre X (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Pre Y (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataYGPU27(para, comm, level);
			 //printf("end exchange Pre Y (level: %d, myID: %d) \n", level, para->getMyID());
			 //printf("start exchange Pre Z (level: %d, myID: %d) \n", level, para->getMyID());
			 exchangePreCollDataZGPU27(para, comm, level);
			 //printf("end exchange Pre Z (level: %d, myID: %d) \n", level, para->getMyID());
			 //////////////////////////////////////////////////////////////////////////
			 //3D domain decomposition convection diffusion
			 if (para->getDiffOn()==true)
			 {
				 exchangePreCollDataADXGPU27(para, comm, level);
				 exchangePreCollDataADYGPU27(para, comm, level);
				 exchangePreCollDataADZGPU27(para, comm, level);
			 }
		 }
		 //////////////////////////////////////////////////////////////////////////////////
		 ////coarse to fine interpolation
			 //ScaleCF_comp_D3Q27F3(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],     para->getParD(level+1)->g6.g[0],
				//			      para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
				//			      para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
				//			      para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
				//			      para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 			   
				//			      para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
				//			      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
				//			      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
				//			      para->getParD(level)->offCF);
    //        getLastCudaError("ScaleCF_comp_D3Q27F3 execution failed");
			////////////////////////////////////////////////////////////////////////
			//ScaleCF_0817_comp_27(para->getParD(level)->d0SP.f[0], para->getParD(level + 1)->d0SP.f[0],
			//				     para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
			//				     para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
			//				     para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
			//				     para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
			//				     para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
			//				     para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
			//				     para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
			//				     para->getParD(level)->offCF);
   //         getLastCudaError("ScaleCF_0817_comp_27 execution failed");
		    ////////////////////////////////////////////////////////////////////////
         ScaleCF_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
								para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
								para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
								para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
								para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								para->getParD(level)->offCF);
         getLastCudaError("ScaleCF27_RhoSq_comp execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCF_RhoSq_3rdMom_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
									//	para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
									//	para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
									//	para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
									//	para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
									//	para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
									//	para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
									//	para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
									//	para->getParD(level)->offCF);
         //getLastCudaError("ScaleCF_RhoSq_3rdMom_comp_27 execution failed");
		 ////////////////////////////////////////////////////////////////////////////////






		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCF27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
         //            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
         //            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP, para->getParD(level+1)->neighborZ_SP,
         //            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
         //            para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
         //            para->getParD(level)->K_CF,           para->getParD(level)->omega,          para->getParD(level+1)->omega, 
         //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny, 
         //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
         //getLastCudaError("ScaleCF27 execution failed");

         //ScaleFC27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
         //            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP, 
         //            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP, para->getParD(level+1)->neighborZ_SP, 
         //            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
         //            para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
         //            para->getParD(level)->K_FC,           para->getParD(level)->omega,          para->getParD(level+1)->omega, 
         //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny, 
         //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
         //getLastCudaError("ScaleFC27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCFEff27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
         //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offCF);
         //getLastCudaError("ScaleCF27 execution failed");

         //ScaleFCEff27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
         //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offFC);
         //getLastCudaError("ScaleFC27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCFLast27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
         //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offCF);
         //getLastCudaError("ScaleCF27 execution failed");

         //ScaleFCLast27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
         //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offFC);
         //getLastCudaError("ScaleFC27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCFpress27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
         //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offCF);
         //getLastCudaError("ScaleCF27 execution failed");

         //ScaleFCpress27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
         //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offFC);
         //getLastCudaError("ScaleFC27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ScaleCF_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
								//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
								//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
								//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								//para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
								//para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								//para->getParD(level)->offCF);
        // getLastCudaError("ScaleCF27_RhoSq_comp execution failed");

        // ScaleFC_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
								//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
								//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
								//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								//para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
								//para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								//para->getParD(level)->offFC);
        // getLastCudaError("ScaleFC27_RhoSq_comp execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // ScaleCF_Fix_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
								//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
								//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
								//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								//para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
								//para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								//para->getParD(level)->offCF);
        // getLastCudaError("ScaleCF27 execution failed");

        // ScaleFC_Fix_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
								//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
								//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
								//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
								//para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
								//para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
								//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
								//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
								//para->getParD(level)->offFC);
        // getLastCudaError("ScaleFC27 execution failed");
         ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
         //ScaleCF_Fix_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
         //               para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offCF);
         //getLastCudaError("ScaleCF27 execution failed");

         //ScaleFC_Fix_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
         //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
         //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
         //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
         //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
         //               para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
         //               para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
         //               para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
         //               para->getParD(level)->offFC);
         //getLastCudaError("ScaleFC27 execution failed");
         //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       //  ScaleCF_NSPress_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
							//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
							//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
							//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
							//para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
							//para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
							//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
							//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
							//para->getParD(level)->offCF);
       //  getLastCudaError("ScaleCF27 execution failed");

       //  ScaleFC_NSPress_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
							//para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
							//para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
							//para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
							//para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
							//para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
							//para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
							//para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
							//para->getParD(level)->offFC);
       //  getLastCudaError("ScaleFC27 execution failed");
		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////
            if (para->getDiffOn())
            {
               if (para->getDiffMod() == 7)
               {
                  //ScaleCFThS7(   para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
                  //               para->getParD(level)->d7.f[0],        para->getParD(level+1)->d7.f[0],                
                  //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                  //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                  //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                  //               para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
                  //               para->getParD(level)->K_CF,           
                  //               para->getParD(level)->vis,            para->getParD(level+1)->diffusivity,   para->getParD(level)->numberofthreads);
                  //getLastCudaError("ScaleCFTh7 execution failed");

                  //ScaleFCThS7(   para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                  //               para->getParD(level)->d7.f[0],        para->getParD(level+1)->d7.f[0],
                  //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
                  //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
                  //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                  //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
                  //               para->getParD(level)->K_FC,
                  //               para->getParD(level)->vis,            para->getParD(level)->diffusivity,     para->getParD(level)->numberofthreads);
                  //getLastCudaError("ScaleFCTh7 execution failed");
                  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  ScaleCFThSMG7( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
                                 para->getParD(level)->d7.f[0],        para->getParD(level+1)->d7.f[0],                
                                 para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                                 para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                                 para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                                 para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
                                 para->getParD(level)->K_CF,           
                                 para->getParD(level)->vis,            para->getParD(level+1)->diffusivity,   para->getParD(level)->numberofthreads,
                                 para->getParD(level)->offCF);
                  getLastCudaError("ScaleCFTh7 execution failed");

                  ScaleFCThSMG7( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                                 para->getParD(level)->d7.f[0],        para->getParD(level+1)->d7.f[0],
                                 para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
                                 para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
                                 para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                                 para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
                                 para->getParD(level)->K_FC,
                                 para->getParD(level)->vis,            para->getParD(level)->diffusivity,     para->getParD(level)->numberofthreads,
                                 para->getParD(level)->offFC);
                  getLastCudaError("ScaleFCTh7 execution failed");
               } 
               else if (para->getDiffMod() == 27)
               {
                  ScaleCFThS27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
                                 para->getParD(level)->d27.f[0],       para->getParD(level+1)->d27.f[0],                
                                 para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                                 para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                                 para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                                 para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF, 
                                 para->getParD(level)->K_CF,           
								 para->getParD(level)->vis,            para->getParD(level+1)->diffusivity,   para->getParD(level)->numberofthreads,
								 para->getParD(level)->offCF);
                  getLastCudaError("ScaleCFTh27 execution failed");

                  ScaleFCThS27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                                 para->getParD(level)->d27.f[0],       para->getParD(level+1)->d27.f[0],
                                 para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
                                 para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
                                 para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                                 para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
                                 para->getParD(level)->K_FC,
								 para->getParD(level)->vis,            para->getParD(level)->diffusivity,     para->getParD(level)->numberofthreads,
								 para->getParD(level)->offFC);
                  getLastCudaError("ScaleFCTh27 execution failed");
               }
            } 
            //////////////////////////////////////////////////////////////////////////

      }
   }
}



