#include "K17CompressibleNavierStokes.h"
#include <logger/Logger.h>
#include "Parameter/Parameter.h"
#include "Parameter/CudaStreamManager.h"
#include "K17CompressibleNavierStokes_Device.cuh"

#include <cuda.h>

template<TurbulenceModel turbulenceModel>
std::shared_ptr< K17CompressibleNavierStokes<turbulenceModel> > K17CompressibleNavierStokes<turbulenceModel>::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<K17CompressibleNavierStokes<turbulenceModel> >(new K17CompressibleNavierStokes<turbulenceModel>(para,level));
}

template<TurbulenceModel turbulenceModel>
void K17CompressibleNavierStokes<turbulenceModel>::run()
{
    K17CompressibleNavierStokes_Device < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads >>>(   para->getParD(level)->omega,
                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                        para->getParD(level)->rho,
                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                        para->getSGSConstant(),
                                                                                                        para->getParD(level)->numberOfNodes,
                                                                                                        level,
                                                                                                        para->getForcesDev(),
                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                        para->getQuadricLimitersDev(),
                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                        para->getParD(level)->taggedFluidNodeIndices[CollisionTemplate::Default],
                                                                                                        para->getParD(level)->numberOfTaggedFluidNodes[CollisionTemplate::Default]);

    getLastCudaError("K17CompressibleNavierStokes_Device execution failed");
}

template<TurbulenceModel turbulenceModel>
void K17CompressibleNavierStokes<turbulenceModel>::runOnIndices( const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex )
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    switch (collisionTemplate)
    {
        case CollisionTemplate::Default:
            K17CompressibleNavierStokes_Device < turbulenceModel, false, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;

        case CollisionTemplate::WriteMacroVars:
            K17CompressibleNavierStokes_Device < turbulenceModel, true, false  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        para->getParD(level)->numberOfNodes,
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
            K17CompressibleNavierStokes_Device < turbulenceModel, true, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>(  para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;

        case CollisionTemplate::ApplyBodyForce:
            K17CompressibleNavierStokes_Device < turbulenceModel, false, true  > <<< cudaGrid.grid, cudaGrid.threads, 0, stream >>>( para->getParD(level)->omega,
                                                                                                                        para->getParD(level)->neighborX, para->getParD(level)->neighborY, para->getParD(level)->neighborZ,
                                                                                                                        para->getParD(level)->distributions.f[0],
                                                                                                                        para->getParD(level)->rho,
                                                                                                                        para->getParD(level)->velocityX, para->getParD(level)->velocityY, para->getParD(level)->velocityZ,
                                                                                                                        para->getParD(level)->turbViscosity,
                                                                                                                        para->getSGSConstant(),
                                                                                                                        para->getParD(level)->numberOfNodes,
                                                                                                                        level,
                                                                                                                        para->getForcesDev(),
                                                                                                                        para->getParD(level)->forceX_SP, para->getParD(level)->forceY_SP, para->getParD(level)->forceZ_SP,
                                                                                                                        para->getQuadricLimitersDev(),
                                                                                                                        para->getParD(level)->isEvenTimestep,
                                                                                                                        indices,
                                                                                                                        size_indices);
            break;
        default:
            throw std::runtime_error("Invalid CollisionTemplate in CumulantK17::runOnIndices()");
            break;
    }

    getLastCudaError("K17CompressibleNavierStokes_Device execution failed");
}

template<TurbulenceModel turbulenceModel>
K17CompressibleNavierStokes<turbulenceModel>::K17CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;

    myPreProcessorTypes.push_back(InitCompSP27);

    

    this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    this->kernelUsesFluidNodeIndices = true;

    VF_LOG_INFO("Using turbulence model: {}", turbulenceModel);
}

template class K17CompressibleNavierStokes<TurbulenceModel::AMD>;
template class K17CompressibleNavierStokes<TurbulenceModel::Smagorinsky>;
template class K17CompressibleNavierStokes<TurbulenceModel::QR>;
template class K17CompressibleNavierStokes<TurbulenceModel::None>;
