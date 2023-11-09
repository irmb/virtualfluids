#include "UpdateGrid27.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#include "Communication/ExchangeData27.h"
#include "Parameter/CudaStreamManager.h"
#include "KernelManager/BCKernelManager.h"
#include "KernelManager/ADKernelManager.h"
#include "KernelManager/GridScalingKernelManager.h"
#include "TurbulenceModels/TurbulenceModelFactory.h"
#include "Kernel/Kernel.h"
#include "Parameter/ParameterRotatingGrid.h"

#include "CollisionStrategy.h"
#include "RefinementStrategy.h"
#include "constants/NumericConstants.h"

#include <parallel/Communicator.h>

void UpdateGrid27::updateGrid(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    SPtr<ParameterRotatingGrid> paraRot = para->getRotatingGridParameter();

    if (level != para->getFine()) {
        updateGrid(level + 1, t);
        if (paraRot == nullptr) updateGrid(level + 1, t);
    }

    //////////////////////////////////////////////////////////////////////////
    
    interactWithProbes(level, t);

    //////////////////////////////////////////////////////////////////////////

    collision(this, para.get(), level, t);

    //////////////////////////////////////////////////////////////////////////

    postCollisionBC(level, t);

    //////////////////////////////////////////////////////////////////////////

    swapBetweenEvenAndOddTimestep(level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getUseWale()) //TODO: make WALE consistent with structure of other turbulence models
        calcMacroscopicQuantities(level);

    calcTurbulentViscosity(level);

    //////////////////////////////////////////////////////////////////////////

    this->preCollisionBC(level, t);

    //////////////////////////////////////////////////////////////////////////
    if (level != para->getFine()) {
        if ( paraRot != nullptr) {
            // rotation
            paraRot->parameterRotHost->gridAngle[0] += paraRot->parameterRotHost->angularVelocity[0];
            paraRot->parameterRotHost->gridAngle[1] += paraRot->parameterRotHost->angularVelocity[1];
            paraRot->parameterRotHost->gridAngle[2] += paraRot->parameterRotHost->angularVelocity[2];
            paraRot->parameterRotDevice->gridAngle[0] += paraRot->parameterRotDevice->angularVelocity[0];
            paraRot->parameterRotDevice->gridAngle[1] += paraRot->parameterRotDevice->angularVelocity[1];
            paraRot->parameterRotDevice->gridAngle[2] += paraRot->parameterRotDevice->angularVelocity[2];
            rotationInterpolation(level);
        } else {
            refinement(this, para.get(), level);
        }
    }

    //////////////////////////////////////////////////////////////////////////
    
    interactWithActuators(level, t);

}

void UpdateGrid27::rotationInterpolation(int level)
{
    // base to nested
    InterpolateStaticToRotating(
        para->getParD(level).get(),
        para->getParD(level+1).get(),
        para->getRotatingGridParameter()->parameterRotDevice.get(),
        &para->getParD(level)->coarseToFine,
        para->getParD(level)->neighborCoarseToFine);

    // nested to base
    InterpolateRotatingToStatic(
        para->getParD(level).get(),
        para->getParD(level+1).get(),
        para->getRotatingGridParameter()->parameterRotDevice.get(),
        &para->getParD(level)->fineToCoarse,
        para->getParD(level)->neighborFineToCoarse);
}

void UpdateGrid27::collisionAllNodes(int level, unsigned int t)
{
    kernels.at(level)->run();

    //////////////////////////////////////////////////////////////////////////

    if (para->getSimulatePorousMedia())
        collisionPorousMedia(level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(level);
}

void UpdateGrid27::collisionUsingIndices(int level, unsigned int t, uint *taggedFluidNodeIndices, uint numberOfTaggedFluidNodes, CollisionTemplate collisionTemplate, CudaStreamIndex stream)
{
    if (taggedFluidNodeIndices != nullptr && numberOfTaggedFluidNodes != 0)
        kernels.at(level)->runOnIndices(taggedFluidNodeIndices, numberOfTaggedFluidNodes, collisionTemplate, stream);
    else
        std::cout << "In collision: fluidNodeIndices or numberOfFluidNodes not defined"
                      << std::endl;

    //////////////////////////////////////////////////////////////////////////
    //! \todo: AD collision and porousMedia should be called separately, not in collisionUsingIndices

    if (para->getSimulatePorousMedia())
        collisionPorousMedia(level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(level);
}

void UpdateGrid27::collisionPorousMedia(int level)
{
    for( std::size_t i = 0; i < pm.size(); i++ )
    {
        KernelPMCumOneCompSP27(para->getParD(level)->numberofthreads,
                               para->getParD(level)->omega,
                               para->getParD(level)->neighborX,
                               para->getParD(level)->neighborY,
                               para->getParD(level)->neighborZ,
                               para->getParD(level)->distributions.f[0],
                               para->getParD(level)->numberOfNodes,
                               level,
                               para->getForcesDev(),
                               pm[i]->getPorosity(),
                               pm[i]->getDarcyLBM(),
                               pm[i]->getForchheimerLBM(),
                               pm[i]->getSizePM(),
                               pm[i]->getHostNodeIDsPM(),
                               para->getParD(level)->isEvenTimestep);
        getLastCudaError("KernelPMCumOneCompSP27 execution failed");
    }
}

void UpdateGrid27::collisionAdvectionDiffusion(int level)
{
    this->adKernelManager->runADcollisionKernel(level);
}

void UpdateGrid27::prepareExchangeMultiGPU(int level, CudaStreamIndex streamIndex)
{
    prepareExchangeCollDataXGPU27AllNodes(para.get(), level, streamIndex);
    prepareExchangeCollDataYGPU27AllNodes(para.get(), level, streamIndex);
    prepareExchangeCollDataZGPU27AllNodes(para.get(), level, streamIndex);
}

void UpdateGrid27::prepareExchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex)
{
    prepareExchangeCollDataXGPU27AfterFtoC(para.get(), level, streamIndex);
    prepareExchangeCollDataYGPU27AfterFtoC(para.get(), level, streamIndex);
    prepareExchangeCollDataZGPU27AfterFtoC(para.get(), level, streamIndex);
}

void UpdateGrid27::exchangeMultiGPU(int level, CudaStreamIndex streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataYGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataZGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AllNodes(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AllNodes(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AllNodes(para.get(), level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }

    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // 1D domain decomposition
    // exchangePostCollDataGPU27(para, comm, level);
}
void UpdateGrid27::exchangeMultiGPU_noStreams_withPrepare(int level, bool useReducedComm)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    if (useReducedComm) {
        // X
        prepareExchangeCollDataXGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataXGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferXGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        // Y
        prepareExchangeCollDataYGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataYGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferYGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        // Z
        prepareExchangeCollDataZGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataZGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferZGPU27AfterFtoC(para.get(), level, CudaStreamIndex::Legacy);
    } else {
        // X
        prepareExchangeCollDataXGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataXGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferXGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        // Y
        prepareExchangeCollDataYGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataYGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferYGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        // Z
        prepareExchangeCollDataZGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
        exchangeCollDataZGPU27AllNodes(para.get(), comm, cudaMemoryManager.get(), level, CudaStreamIndex::Legacy);
        scatterNodesFromRecvBufferZGPU27AllNodes(para.get(), level, CudaStreamIndex::Legacy);
    }

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }
}
void UpdateGrid27::exchangeMultiGPUAfterFtoC(int level, CudaStreamIndex streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataYGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);
    exchangeCollDataZGPU27AfterFtoC(para.get(), comm, cudaMemoryManager.get(), level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AfterFtoC(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AfterFtoC(para.get(), level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AfterFtoC(para.get(), level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADYGPU27(para.get(), comm, cudaMemoryManager.get(), level);
        exchangePostCollDataADZGPU27(para.get(), comm, cudaMemoryManager.get(), level);
    }
}

void UpdateGrid27::postCollisionBC(int level, uint t)
{
    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    // V E L O C I T Y (I N F L O W)
    this->bcKernelManager->runVelocityBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // N O - S L I P
    this->bcKernelManager->runNoSlipBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // S L I P
    this->bcKernelManager->runSlipBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // S T R E S S (wall model)
    this->bcKernelManager->runStressWallModelKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    this->bcKernelManager->runGeoBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    this->bcKernelManager->runOutflowBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    this->bcKernelManager->runPressureBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // P R E C U R S O R
    this->bcKernelManager->runPrecursorBCKernelPost(level, t, cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    if (para->getDiffOn())
    {
        this->adKernelManager->runADgeometryBCKernel(level);
        this->adKernelManager->runADveloBCKernel(level);
        this->adKernelManager->runADslipBCKernel(level);
        this->adKernelManager->runADpressureBCKernel(level);
    }
}

void UpdateGrid27::swapBetweenEvenAndOddTimestep(int level)
{
    if (para->getParD(level)->isEvenTimestep==true)  para->getParD(level)->isEvenTimestep=false;
    else                                        para->getParD(level)->isEvenTimestep=true;
}

void UpdateGrid27::calcMacroscopicQuantities(int level)
{
    CalcMacCompSP27(para->getParD(level)->velocityX,
                    para->getParD(level)->velocityY,
                    para->getParD(level)->velocityZ,
                    para->getParD(level)->rho,
                    para->getParD(level)->pressure,
                    para->getParD(level)->typeOfGridNode,
                    para->getParD(level)->neighborX,
                    para->getParD(level)->neighborY,
                    para->getParD(level)->neighborZ,
                    para->getParD(level)->numberOfNodes,
                    para->getParD(level)->numberofthreads,
                    para->getParD(level)->distributions.f[0],
                    para->getParD(level)->isEvenTimestep);
    getLastCudaError("CalcMacSP27 execution failed");
}

void UpdateGrid27::preCollisionBC(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    // V E L O C I T Y (I N F L O W)
    this->bcKernelManager->runVelocityBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    this->bcKernelManager->runGeoBCKernelPre(level, t, cudaMemoryManager.get());

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    this->bcKernelManager->runPressureBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    this->bcKernelManager->runOutflowBCKernelPre(level);

    //////////////////////////////////////////////////////////////////////////////////
    ////only for a round off error test
    //para->cudaCopyTestREtoHost(0,para->getParH(0)->pressureBC.numberOfBCnodes);
    //printRE(para, t);
    //////////////////////////////////////////////////////////////////////////////////
}

void UpdateGrid27::fineToCoarse(int level, InterpolationCells* fineToCoarse, ICellNeigh &neighborFineToCoarse, CudaStreamIndex streamIndex)
{
    gridScalingKernelManager->runFineToCoarseKernelLB(level, fineToCoarse, neighborFineToCoarse, streamIndex);

    if (para->getDiffOn()) {
        if (para->getStreamManager()->streamIsRegistered(streamIndex)) {
            printf("fineToCoarse Advection Diffusion not implemented"); // TODO
            return;
        }
        gridScalingKernelManager->runFineToCoarseKernelAD(level);
    }
}

void UpdateGrid27::coarseToFine(int level, InterpolationCells* coarseToFine, ICellNeigh &neighborCoarseToFine, CudaStreamIndex streamIndex)
{
    this->gridScalingKernelManager->runCoarseToFineKernelLB(level, coarseToFine, neighborCoarseToFine, streamIndex);

    if (para->getDiffOn())
    {
        if(para->getStreamManager()->streamIsRegistered(streamIndex)){
            printf("CoarseToFineWithStream Advection Diffusion not implemented"); // TODO
            return;
        }
        this->gridScalingKernelManager->runCoarseToFineKernelAD(level);
    }
}

void UpdateGrid27::interactWithActuators(int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> actuator: para->getActuators() )
    {
        actuator->interact(para.get(), cudaMemoryManager.get(), level, t);
    }
}

void  UpdateGrid27::interactWithProbes(int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> probe: para->getProbes() )
    {
        probe->interact(para.get(), cudaMemoryManager.get(), level, t);
    }
}

void  UpdateGrid27::calcTurbulentViscosity(int level)
{
    this->tmFactory->runTurbulenceModelKernel(level);
}

void UpdateGrid27::exchangeData(int level)
{
    exchangeMultiGPU_noStreams_withPrepare(level, false);
}

UpdateGrid27::UpdateGrid27(SPtr<Parameter> para, vf::parallel::Communicator &comm, SPtr<CudaMemoryManager> cudaMemoryManager,
                           std::vector<std::shared_ptr<PorousMedia>> &pm, std::vector<SPtr<Kernel>> &kernels,
                           BoundaryConditionFactory *bcFactory, SPtr<TurbulenceModelFactory> tmFactory,
                           GridScalingFactory *scalingFactory)
    : para(para), comm(comm), cudaMemoryManager(cudaMemoryManager), pm(pm), kernels(kernels), tmFactory(tmFactory)
{
    this->collision = getFunctionForCollisionAndExchange(para->getUseStreams(), para->getNumprocs(), para->getKernelNeedsFluidNodeIndicesToRun());
    this->refinement = getFunctionForRefinementAndExchange(para->getUseStreams(), para->getNumprocs(), para->getMaxLevel(), para->useReducedCommunicationAfterFtoC);

    this->bcKernelManager = std::make_shared<BCKernelManager>(para, bcFactory);
    this->adKernelManager = std::make_shared<ADKernelManager>(para);
    this->gridScalingKernelManager = std::make_shared<GridScalingKernelManager>(para, scalingFactory);

    SPtr<ParameterRotatingGrid> paraRot = para->getRotatingGridParameter();

    uint level = 0;
    real finalRotation = (real) (vf::basics::constant::cPi);
    VF_LOG_INFO("intended rotation: {} around center point ({}, {}, {})", finalRotation,
                paraRot->parameterRotHost->centerPoint[0], paraRot->parameterRotHost->centerPoint[1],
                paraRot->parameterRotHost->centerPoint[2]);
    real steps = 800;
    paraRot->parameterRotHost->angularVelocity[0] = finalRotation / steps;
    paraRot->parameterRotDevice->angularVelocity[0] = paraRot->parameterRotHost->angularVelocity[0];
    VF_LOG_INFO("angular velocity: {}", paraRot->parameterRotHost->angularVelocity[0]);
    for (int i = 0; i < steps; i++) {
        paraRot->parameterRotHost->gridAngle[0] += paraRot->parameterRotHost->angularVelocity[0];
        paraRot->parameterRotHost->gridAngle[1] += paraRot->parameterRotHost->angularVelocity[1];
        paraRot->parameterRotHost->gridAngle[2] += paraRot->parameterRotHost->angularVelocity[2];
        paraRot->parameterRotDevice->gridAngle[0] += paraRot->parameterRotDevice->angularVelocity[0];
        paraRot->parameterRotDevice->gridAngle[1] += paraRot->parameterRotDevice->angularVelocity[1];
        paraRot->parameterRotDevice->gridAngle[2] += paraRot->parameterRotDevice->angularVelocity[2];

        // base to nested
        TraverseStaticToRotating(
            para->getParD(level).get(),
            para->getParD(level+1).get(),
            para->getRotatingGridParameter()->parameterRotDevice.get(),
            &para->getParD(level)->coarseToFine,
            para->getParD(level)->neighborCoarseToFine);

        // nested to base
        TraverseRotatingToStatic(
            para->getParD(level).get(),
            para->getParD(level+1).get(),
            para->getRotatingGridParameter()->parameterRotDevice.get(),
            &para->getParD(level)->fineToCoarse,
            para->getParD(level)->neighborFineToCoarse);
    }
    paraRot->parameterRotHost->angularVelocity[0] = 0.0;
    paraRot->parameterRotHost->angularVelocity[1] = 0.0;
    paraRot->parameterRotHost->angularVelocity[2] = 0.0;
    paraRot->parameterRotDevice->angularVelocity[0] = 0.0;
    paraRot->parameterRotDevice->angularVelocity[1] = 0.0;
    paraRot->parameterRotDevice->angularVelocity[2] = 0.0;

    VF_LOG_INFO("Angle x {}", paraRot->parameterRotHost->gridAngle[0]);
    VF_LOG_INFO("Angle y {}", paraRot->parameterRotHost->gridAngle[1]);
    VF_LOG_INFO("Angle z {}", paraRot->parameterRotHost->gridAngle[2]);
    VF_LOG_INFO("Angular velocity x {}", paraRot->parameterRotHost->angularVelocity[0]);
    VF_LOG_INFO("Angular velocity y {}", paraRot->parameterRotHost->angularVelocity[1]);
    VF_LOG_INFO("Angular velocity z {}", paraRot->parameterRotHost->angularVelocity[2]);
}
