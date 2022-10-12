#include "CumulantK17Almighty.h"
#include <logger/Logger.h>
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "CumulantK17Almighty_Device.cuh"

#include <cuda.h>

template<TurbulenceModel turbulenceModel> 
std::shared_ptr< CumulantK17Almighty<turbulenceModel> > CumulantK17Almighty<turbulenceModel>::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
	return std::shared_ptr<CumulantK17Almighty<turbulenceModel> >(new CumulantK17Almighty<turbulenceModel>(para,level));
}

template<TurbulenceModel turbulenceModel>
void CumulantK17Almighty<turbulenceModel>::run()
{
	LB_Kernel_CumulantK17Almighty < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads >>>(   para->getParD(level)->omega, 	
																												para->getParD(level)->typeOfGridNode, 										
																												para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,	
																												para->getParD(level)->distributions.f[0],	
																												para->getParD(level)->rho,		
																												para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,	
																												para->getParD(level)->turbViscosity,
																												para->getSGSConstant(),
																												(unsigned long)para->getParD(level)->numberOfNodes,	
																												level,				
																												para->getIsBodyForce(),				
																												para->getForcesDev(),				
																												para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
																												para->getQuadricLimitersDev(),			
																												para->getParD(level)->isEvenTimestep,
																												para->getParD(level)->fluidNodeIndices,
        																										para->getParD(level)->numberOfFluidNodes);

	getLastCudaError("LB_Kernel_CumulantK17Almighty execution failed");
}

template<TurbulenceModel turbulenceModel>
void CumulantK17Almighty<turbulenceModel>::runOnIndices( const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, int streamIndex )
{
	cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager()->getStream(streamIndex);

	switch (collisionTemplate)
	{
		case CollisionTemplate::Default:
			LB_Kernel_CumulantK17Almighty < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(   	para->getParD(level)->omega, 	
																																	para->getParD(level)->typeOfGridNode, 										
																																	para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,	
																																	para->getParD(level)->distributions.f[0],	
																																	para->getParD(level)->rho,		
																																	para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,	
																																	para->getParD(level)->turbViscosity,
																																	para->getSGSConstant(),
																																	(unsigned long)para->getParD(level)->numberOfNodes,	
																																	level,				
																																	para->getIsBodyForce(),				
																																	para->getForcesDev(),				
																																	para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
																																	para->getQuadricLimitersDev(),			
																																	para->getParD(level)->isEvenTimestep,
																																	indices,
																																	size_indices);
			break;
		
		case CollisionTemplate::WriteMacroVars:
			LB_Kernel_CumulantK17Almighty < turbulenceModel, true, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega, 	
																																para->getParD(level)->typeOfGridNode, 										
																																para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,	
																																para->getParD(level)->distributions.f[0],	
																																para->getParD(level)->rho,		
																																para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,	
																																para->getParD(level)->turbViscosity,
																																para->getSGSConstant(),
																																(unsigned long)para->getParD(level)->numberOfNodes,	
																																level,				
																																para->getIsBodyForce(),				
																																para->getForcesDev(),				
																																para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
																																para->getQuadricLimitersDev(),			
																																para->getParD(level)->isEvenTimestep,
																																indices,
																																size_indices);
			break;
		
		case CollisionTemplate::AllFeatures:
			LB_Kernel_CumulantK17Almighty < turbulenceModel, true, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(  para->getParD(level)->omega, 	
																																para->getParD(level)->typeOfGridNode, 										
																																para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,	
																																para->getParD(level)->distributions.f[0],	
																																para->getParD(level)->rho,		
																																para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,	
																																para->getParD(level)->turbViscosity,
																																para->getSGSConstant(),
																																(unsigned long)para->getParD(level)->numberOfNodes,	
																																level,				
																																para->getIsBodyForce(),				
																																para->getForcesDev(),				
																																para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
																																para->getQuadricLimitersDev(),			
																																para->getParD(level)->isEvenTimestep,
																																indices,
																																size_indices);
			break;
		case CollisionTemplate::ApplyBodyForce:
			LB_Kernel_CumulantK17Almighty < turbulenceModel, false, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega, 	
																																para->getParD(level)->typeOfGridNode, 										
																																para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,	
																																para->getParD(level)->distributions.f[0],	
																																para->getParD(level)->rho,		
																																para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,	
																																para->getParD(level)->turbViscosity,
																																para->getSGSConstant(),
																																(unsigned long)para->getParD(level)->numberOfNodes,	
																																level,				
																																para->getIsBodyForce(),				
																																para->getForcesDev(),				
																																para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
																																para->getQuadricLimitersDev(),			
																																para->getParD(level)->isEvenTimestep,
																																indices,
																																size_indices);
			break;
		default:
			throw std::runtime_error("Invalid CollisionTemplate in CumulantK17Almighty::runOnIndices()");
			break;
	}

	getLastCudaError("LB_Kernel_CumulantK17Almighty execution failed");
}

template<TurbulenceModel turbulenceModel>
CumulantK17Almighty<turbulenceModel>::CumulantK17Almighty(std::shared_ptr<Parameter> para, int level)
{
	this->para = para;
	this->level = level;

	myPreProcessorTypes.push_back(InitCompSP27);

	myKernelGroup = BasicKernel;

	this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
	this->kernelUsesFluidNodeIndices = true;
	
	VF_LOG_INFO("Using turbulence model: {}", turbulenceModel);
}

template class CumulantK17Almighty<TurbulenceModel::AMD>;
template class CumulantK17Almighty<TurbulenceModel::Smagorinsky>;
template class CumulantK17Almighty<TurbulenceModel::QR>;
template class CumulantK17Almighty<TurbulenceModel::None>;
