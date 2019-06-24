#include "Init/InitLattice.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Parameter/Parameter.h"
#include "GPU/GPU_Interface.h"
#include "Temperature/FindTemperature.h"
#include "PreProcessor/PreProcessor.h"
#include "GPU/CudaMemoryManager.h"

////////////////////////////////////////////////////////////////////////////////
void initLattice(SPtr<Parameter> para, SPtr<PreProcessor> preProcessor, SPtr<CudaMemoryManager> cudaManager)
{
    for (int lev=para->getFine(); lev >= para->getCoarse(); lev--)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		preProcessor->init(para, lev);

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //CalcMacSP27(para->getParD(lev)->vx_SP,       
        //            para->getParD(lev)->vy_SP,        
        //            para->getParD(lev)->vz_SP,        
        //            para->getParD(lev)->rho_SP, 
        //            para->getParD(lev)->press_SP, 
        //            para->getParD(lev)->geoSP,       
        //            para->getParD(lev)->neighborX_SP, 
        //            para->getParD(lev)->neighborY_SP, 
        //            para->getParD(lev)->neighborZ_SP,
        //            para->getParD(lev)->size_Mat_SP, 
        //            para->getParD(lev)->numberofthreads,       
        //            para->getParD(lev)->d0SP.f[0],    
        //            para->getParD(lev)->evenOrOdd);
        //getLastCudaError("Kernel CalcMacSP27 execution failed"); 
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
        getLastCudaError("Kernel execution failed"); 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (para->getCalcMedian())
		{
			unsigned int tdiff = 1;
			CalcMacMedSP27( para->getParD(lev)->vx_SP_Med,       
							para->getParD(lev)->vy_SP_Med,        
							para->getParD(lev)->vz_SP_Med,        
							para->getParD(lev)->rho_SP_Med, 
							para->getParD(lev)->press_SP_Med, 
							para->getParD(lev)->geoSP,       
							para->getParD(lev)->neighborX_SP, 
							para->getParD(lev)->neighborY_SP, 
							para->getParD(lev)->neighborZ_SP,
							tdiff,
							para->getParD(lev)->size_Mat_SP, 
							para->getParD(lev)->numberofthreads,       
							para->getParD(lev)->evenOrOdd);
			getLastCudaError("CalcMacMedSP27 execution failed"); 
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// advection - diffusion stuff
		if (para->getDiffOn()==true){
			//malloc
			//printf("vor cudaAllocConc\n");
			cudaManager->cudaAllocConc(lev);
			//define init conc/temp
			//printf("vor Schleife\n");
			for (unsigned int i = 0; i < para->getParH(lev)->size_Mat_SP; i++)
			{
				para->getParH(lev)->Conc[i] = para->getTemperatureInit();
			}
			//malloc and init fs
			//printf("vor initTemperatur\n");
			initTemperatur(para.get(), cudaManager.get(), lev);
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}
