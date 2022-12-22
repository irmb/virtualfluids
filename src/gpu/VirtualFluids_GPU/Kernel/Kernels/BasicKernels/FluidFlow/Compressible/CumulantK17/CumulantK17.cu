#include "CumulantK17.h"
#include <logger/Logger.h>
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "CumulantK17_Device.cuh"

#include <cuda.h>

template<TurbulenceModel turbulenceModel>
std::shared_ptr< CumulantK17<turbulenceModel> > CumulantK17<turbulenceModel>::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<CumulantK17<turbulenceModel> >(new CumulantK17<turbulenceModel>(para,level));
}

template<TurbulenceModel turbulenceModel>
void CumulantK17<turbulenceModel>::run()
{
    LB_Kernel_CumulantK17 < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads >>>(   para->getParD(level)->omega,
                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                        para->getParD(level)->rho,
                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                        para->getSGSConstant(),
                                                                                                        (unsigned long)para->getParD(level)->numberOfNodes,
                                                                                                        level,
                                                                                                        para->getForcesDev(),
                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                        para->getQuadricLimitersDev(),
                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                        para->getParD(level)->taggedFluidNodeIndices[CollisionTemplate::Default],
                                                                                                        para->getParD(level)->numberOfTaggedFluidNodes[CollisionTemplate::Default]);

    getLastCudaError("LB_Kernel_CumulantK17 execution failed");
}

template<TurbulenceModel turbulenceModel>
void CumulantK17<turbulenceModel>::runOnIndices( const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex )
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    switch (collisionTemplate)
    {
        case CollisionTemplate::Default:
            LB_Kernel_CumulantK17 < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        (unsigned long)para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;

        case CollisionTemplate::WriteMacroVars:
            LB_Kernel_CumulantK17 < turbulenceModel, true, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        (unsigned long)para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;

        case CollisionTemplate::SubDomainBorder:
        case CollisionTemplate::AllFeatures:
            LB_Kernel_CumulantK17 < turbulenceModel, true, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(  para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        (unsigned long)para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;	case CollisionTemplate::ApplyBodyForce:
            LB_Kernel_CumulantK17 < turbulenceModel, false, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        (unsigned long)para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;	default:
            throw std::runtime_error("Invalid CollisionTemplate in CumulantK17::runOnIndices()");
            break;
    }

    getLastCudaError("LB_Kernel_CumulantK17 execution failed");
}

template<TurbulenceModel turbulenceModel>
CumulantK17<turbulenceModel>::CumulantK17(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);

    myKernelGroup = BasicKernel;

    this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    this->kernelUsesFluidNodeIndices = true;

    VF_LOG_INFO("Using turbulence model: {}", turbulenceModel);
}

template class CumulantK17<TurbulenceModel::AMD>;
template class CumulantK17<TurbulenceModel::Smagorinsky>;
template class CumulantK17<TurbulenceModel::QR>;
template class CumulantK17<TurbulenceModel::None>;
