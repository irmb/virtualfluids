#include "Init/InitLattice.h"

#include "LBM/LB.h"
#include "GPU/GPU_Interface.h"

#include <cuda_runtime.h>

#include "Parameter/Parameter.h"

void initLattice(std::shared_ptr<Parameter> para)
{
  para->getParD(0)->evenOrOdd = false;

  InitSP27(para->getParD(0)->numberofthreads, 
           para->getParD(0)->neighborX_SP, 
           para->getParD(0)->neighborY_SP, 
           para->getParD(0)->neighborZ_SP, 
           para->getParD(0)->geoSP,       
           para->getParD(0)->rho_SP, 
           para->getParD(0)->vx_SP, 
           para->getParD(0)->vy_SP, 
           para->getParD(0)->vz_SP, 
           para->getParD(0)->size_Mat_SP, 
           para->getParD(0)->d0SP.f[0],
           para->getParD(0)->evenOrOdd); 
  getLastCudaError("Kernel execution failed"); 

  para->getParD(0)->evenOrOdd = true;

  InitSP27(para->getParD(0)->numberofthreads, 
           para->getParD(0)->neighborX_SP, 
           para->getParD(0)->neighborY_SP, 
           para->getParD(0)->neighborZ_SP, 
           para->getParD(0)->geoSP,       
           para->getParD(0)->rho_SP, 
           para->getParD(0)->vx_SP, 
           para->getParD(0)->vy_SP, 
           para->getParD(0)->vz_SP, 
           para->getParD(0)->size_Mat_SP, 
           para->getParD(0)->d0SP.f[0],
           para->getParD(0)->evenOrOdd); 
  getLastCudaError("Kernel execution failed"); 

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
  getLastCudaError("Kernel CalcMacCompSP27 execution failed"); 
}
