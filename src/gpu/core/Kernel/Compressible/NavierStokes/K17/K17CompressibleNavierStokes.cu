#include "K17CompressibleNavierStokes.h"

#include <optional>
#include <stdexcept>

#include <cuda.h>

#include <logger/Logger.h>

#include <lbm/collision/TurbulentViscosity.h>
#include <lbm/collision/K17CompressibleNavierStokes.h>

#include "Kernel/KernelImp.h"
#include "Cuda/CudaStreamManager.h"
#include "Parameter/Parameter.h"
#include "LBM/GPUHelperFunctions/KernelUtilities.h"
#include "LBM/GPUHelperFunctions/RunCollision.cuh"

template <vf::lbm::TurbulenceModel turbulenceModel>
std::shared_ptr<K17CompressibleNavierStokes<turbulenceModel>>
K17CompressibleNavierStokes<turbulenceModel>::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<K17CompressibleNavierStokes<turbulenceModel>>(new K17CompressibleNavierStokes(para, level));
}

vf::gpu::GPUCollisionParameter getCollisionParameter(const std::shared_ptr<Parameter>& para, int level, const uint* indices,
                                               uint sizeIndices)
{
    vf::gpu::GPUCollisionParameter collisionParameter { para->getParD(level)->omega,
                                                        para->getParD(level)->neighborX,
                                                        para->getParD(level)->neighborY,
                                                        para->getParD(level)->neighborZ,
                                                        para->getParD(level)->distributions.f[0],
                                                        para->getParD(level)->rho,
                                                        para->getParD(level)->velocityX,
                                                        para->getParD(level)->velocityY,
                                                        para->getParD(level)->velocityZ,
                                                        para->getParD(level)->turbViscosity,
                                                        para->getSGSConstant(),
                                                        (int)para->getParD(level)->numberOfNodes,
                                                        vf::gpu::getForceFactor(level),
                                                        para->getForcesDev(),
                                                        para->getParD(level)->forceX_SP,
                                                        para->getParD(level)->forceY_SP,
                                                        para->getParD(level)->forceZ_SP,
                                                        para->getQuadricLimitersDev(),
                                                        para->getParD(level)->isEvenTimestep,
                                                        indices,
                                                        sizeIndices };

    return collisionParameter;
}

template <vf::lbm::TurbulenceModel turbulenceModel>
void K17CompressibleNavierStokes<turbulenceModel>::run()
{
    vf::gpu::GPUCollisionParameter kernelParameter =
        getCollisionParameter(para, level, para->getParD(level)->taggedFluidNodeIndices[CollisionTemplate::Default],
                           para->getParD(level)->numberOfTaggedFluidNodes[CollisionTemplate::Default]);

    auto collision = [] __device__(vf::lbm::CollisionParameter & parameter, vf::lbm::MacroscopicValues & macroscopicValues,
                                   vf::lbm::TurbulentViscosity & turbulentViscosity) {
        return vf::lbm::runK17CompressibleNavierStokes<turbulenceModel>(parameter, macroscopicValues, turbulentViscosity);
    };

    vf::gpu::runCollision<decltype(collision), turbulenceModel, false, false><<<cudaGrid.grid, cudaGrid.threads>>>(collision, kernelParameter);

    getLastCudaError("K17CompressibleNavierStokesUnified execution failed");
}

template <vf::lbm::TurbulenceModel turbulenceModel>
void K17CompressibleNavierStokes<turbulenceModel>::runOnIndices(const unsigned int* indices, unsigned int size_indices,
                                                                CollisionTemplate collisionTemplate,
                                                                CudaStreamIndex streamIndex)
{
    cudaStream_t stream = para->getStreamManager()->getStream(streamIndex);

    auto collision = [] __device__(vf::lbm::CollisionParameter& parameter, vf::lbm::MacroscopicValues& macroscopicValues, vf::lbm::TurbulentViscosity& turbulentViscosity) {
        return vf::lbm::runK17CompressibleNavierStokes<turbulenceModel>(parameter, macroscopicValues, turbulentViscosity);
    };

    switch (collisionTemplate) {
        case CollisionTemplate::Default: {
            vf::gpu::GPUCollisionParameter kernelParameter = getCollisionParameter(para, level, indices, size_indices);
            vf::gpu::runCollision<decltype(collision), turbulenceModel, false, false><<<cudaGrid.grid, cudaGrid.threads, 0, stream>>>(collision, kernelParameter);

            break;
        }
        case CollisionTemplate::WriteMacroVars: {
            vf::gpu::GPUCollisionParameter kernelParameter = getCollisionParameter(para, level, indices, size_indices);
            vf::gpu::runCollision<decltype(collision), turbulenceModel, true, false><<<cudaGrid.grid, cudaGrid.threads, 0, stream>>>(collision, kernelParameter);

            break;
        }
        case CollisionTemplate::SubDomainBorder:
        case CollisionTemplate::AllFeatures: {
            vf::gpu::GPUCollisionParameter kernelParameter = getCollisionParameter(para, level, indices, size_indices);
            vf::gpu::runCollision<decltype(collision), turbulenceModel, true, true><<<cudaGrid.grid, cudaGrid.threads, 0, stream>>>(collision, kernelParameter);

            break;
        }
        case CollisionTemplate::ApplyBodyForce: {
            vf::gpu::GPUCollisionParameter kernelParameter = getCollisionParameter(para, level, indices, size_indices);
            vf::gpu::runCollision<decltype(collision), turbulenceModel, false, true><<<cudaGrid.grid, cudaGrid.threads, 0, stream>>>(collision, kernelParameter);

            break;
        }
        default:
            throw std::runtime_error("Invalid CollisionTemplate in CumulantK17::runOnIndices()");
            break;
    }

    getLastCudaError("K17CompressibleNavierStokes_Device execution failed");
}

template <vf::lbm::TurbulenceModel turbulenceModel>
K17CompressibleNavierStokes<turbulenceModel>::K17CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level) : KernelImp(para, level)
{
    myPreProcessorTypes.push_back(InitNavierStokesCompressible);

    this->cudaGrid = vf::cuda::CudaGrid(para->getParD(level)->numberofthreads, para->getParD(level)->numberOfNodes);
    this->kernelUsesFluidNodeIndices = true;

    VF_LOG_INFO("Using turbulence model: {}", turbulenceModel);
}

template class K17CompressibleNavierStokes<vf::lbm::TurbulenceModel::AMD>;
template class K17CompressibleNavierStokes<vf::lbm::TurbulenceModel::Smagorinsky>;
template class K17CompressibleNavierStokes<vf::lbm::TurbulenceModel::QR>;
template class K17CompressibleNavierStokes<vf::lbm::TurbulenceModel::None>;
