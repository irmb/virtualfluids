#include "UpdateGrid27.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "Calculation/DragLift.h"
#include "Calculation/Cp.h"
//#include "Utilities/StringUtil.hpp"
//#include "Output/UnstructuredGridWriter.hpp"
#include "Communication/ExchangeData27.h"
#include "Kernel/Kernel.h"
#include "Parameter/CudaStreamManager.h"
#include "GPU/TurbulentViscosity.h"
#include "GPU/CudaKernelManager.h"

void UpdateGrid27::updateGrid(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////

    if (level != para->getFine()) {
        updateGrid(level + 1, t);
        updateGrid(level + 1, t);
    }

    //////////////////////////////////////////////////////////////////////////

    (this->*collisionAndExchange)(level, t);

    //////////////////////////////////////////////////////////////////////////

    this->postCollisionBC(level, t);

    //////////////////////////////////////////////////////////////////////////

    swapBetweenEvenAndOddTimestep(para.get(), level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getUseWale())
        calcMacroscopicQuantities(para.get(), level);

    if (para->getUseTurbulentViscosity())
        calcTurbulentViscosity(para.get(), level);

    //////////////////////////////////////////////////////////////////////////

    preCollisionBC(para.get(), cudaManager.get(), level, t);

    //////////////////////////////////////////////////////////////////////////
    if( level != para->getFine() )
    {
        (this->*refinementAndExchange)(level);
    }

    interactWithActuators(para.get(), cudaManager.get(), level, t);

    interactWithProbes(para.get(), cudaManager.get(), level, t);
}

void UpdateGrid27::refinementAndExchange_noRefinementAndExchange(int level) {}

void UpdateGrid27::refinementAndExchange_streams_onlyExchangeInterface(int level)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex   = para->getStreamManager()->getBulkStreamIndex();

    // fine to coarse border
    fineToCoarseWithStream(para.get(), level, para->getParD(level)->intFCBorder.ICellFCC,
                           para->getParD(level)->intFCBorder.ICellFCF, para->getParD(level)->intFCBorder.kFC,
                           borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    prepareExchangeMultiGPUAfterFtoC(para.get(), level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernels (f to c and c to f)
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    fineToCoarseWithStream(para.get(), level, para->getParD(level)->intFCBulk.ICellFCC,
                           para->getParD(level)->intFCBulk.ICellFCF, para->getParD(level)->intFCBulk.kFC,
                           bulkStreamIndex);
    coarseToFineWithStream(para.get(), level, para->getParD(level)->intCFBulk.ICellCFC,
                           para->getParD(level)->intCFBulk.ICellCFF, para->getParD(level)->intCFBulk.kCF, para->getParD(level)->offCFBulk,
                           bulkStreamIndex);

    // exchange
    exchangeMultiGPUAfterFtoC(para.get(), comm, cudaManager.get(), level, borderStreamIndex);

    // coarse to fine border
    coarseToFineWithStream(para.get(), level, para->getParD(level)->intCFBorder.ICellCFC,
                           para->getParD(level)->intCFBorder.ICellCFF, para->getParD(level)->intCFBorder.kCF, para->getParD(level)->offCF,
                           borderStreamIndex);
    cudaDeviceSynchronize();
}

void UpdateGrid27::refinementAndExchange_streams_completeExchange(int level)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex   = para->getStreamManager()->getBulkStreamIndex();

    // fine to coarse border
    fineToCoarseWithStream(para.get(), level, para->getParD(level)->intFCBorder.ICellFCC,
                           para->getParD(level)->intFCBorder.ICellFCF, para->getParD(level)->intFCBorder.kFC,
                           borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    prepareExchangeMultiGPU(para.get(), level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernels (f to c and c to f)
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    fineToCoarseWithStream(para.get(), level, para->getParD(level)->intFCBulk.ICellFCC,
                           para->getParD(level)->intFCBulk.ICellFCF, para->getParD(level)->intFCBulk.kFC,
                           bulkStreamIndex);
    coarseToFineWithStream(para.get(), level, para->getParD(level)->intCFBulk.ICellCFC,
                           para->getParD(level)->intCFBulk.ICellCFF, para->getParD(level)->intCFBulk.kCF, para->getParD(level)->offCFBulk,
                           bulkStreamIndex);

    // exchange
    exchangeMultiGPU(para.get(), comm, cudaManager.get(), level, borderStreamIndex);

    // coarse to fine border
    coarseToFineWithStream(para.get(), level, para->getParD(level)->intCFBorder.ICellCFC,
                           para->getParD(level)->intCFBorder.ICellCFF, para->getParD(level)->intCFBorder.kCF, para->getParD(level)->offCF,
                           borderStreamIndex);
    cudaDeviceSynchronize();
}

void UpdateGrid27::refinementAndExchange_noStreams_onlyExchangeInterface(int level)
{
    fineToCoarse(para.get(), level);

    exchangeMultiGPU_noStreams_withPrepare(para.get(), comm, cudaManager.get(), level, true);

    coarseToFine(para.get(), level);
}

void UpdateGrid27::refinementAndExchange_noStreams_completeExchange(int level)
{
    fineToCoarse(para.get(), level);

    exchangeMultiGPU_noStreams_withPrepare(para.get(), comm, cudaManager.get(), level, false);

    coarseToFine(para.get(), level);
}

void UpdateGrid27::refinementAndExchange_noExchange(int level)
{
    fineToCoarse(para.get(), level);
    coarseToFine(para.get(), level);
}

void UpdateGrid27::collisionAndExchange_noStreams_indexKernel(int level, unsigned int t)
{
    collisionUsingIndex(para.get(), pm, level, t, kernels, para->getParD(level)->fluidNodeIndices,
                            para->getParD(level)->numberOfFluidNodes, -1);
    exchangeMultiGPU_noStreams_withPrepare(para.get(), comm, cudaManager.get(), level, false);
}

void UpdateGrid27::collisionAndExchange_noStreams_oldKernel(int level, unsigned int t)
{
    collision(para.get(), pm, level, t, kernels);
    exchangeMultiGPU_noStreams_withPrepare(para.get(), comm, cudaManager.get(), level, false);
}

void UpdateGrid27::collisionAndExchange_streams(int level, unsigned int t)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex   = para->getStreamManager()->getBulkStreamIndex();
    // launch border kernel
    collisionUsingIndex(para.get(), pm, level, t, kernels, para->getParD(level)->fluidNodeIndicesBorder,
                        para->getParD(level)->numberOffluidNodesBorder, borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    prepareExchangeMultiGPU(para.get(), level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernel
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    collisionUsingIndex(para.get(), pm, level, t, kernels, para->getParD(level)->fluidNodeIndices,
                        para->getParD(level)->numberOfFluidNodes, bulkStreamIndex);

    exchangeMultiGPU(para.get(), comm, cudaManager.get(), level, borderStreamIndex);
}

void collision(Parameter* para, std::vector<std::shared_ptr<PorousMedia>>& pm, int level, unsigned int t, std::vector < SPtr< Kernel>>& kernels)
{
    kernels.at(level)->run();

    //////////////////////////////////////////////////////////////////////////

    if (para->getSimulatePorousMedia())
        collisionPorousMedia(para, pm, level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(para, level);
}

void collisionUsingIndex(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,
                         std::vector<SPtr<Kernel>> &kernels, uint *fluidNodeIndices, uint numberOfFluidNodes, int stream)
{
    if (fluidNodeIndices != nullptr && numberOfFluidNodes != 0)
        kernels.at(level)->runOnIndices(fluidNodeIndices, numberOfFluidNodes, stream);
    else
        std::cout << "In collision: fluidNodeIndices or numberOfFluidNodes not definded"
                      << std::endl;

    //////////////////////////////////////////////////////////////////////////

    if (para->getSimulatePorousMedia())
        collisionPorousMedia(para, pm, level);

    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
        collisionAdvectionDiffusion(para, level);
}

void collisionPorousMedia(Parameter* para, std::vector<std::shared_ptr<PorousMedia>>& pm, int level)
{
    for( std::size_t i = 0; i < pm.size(); i++ )
    {
        KernelPMCumOneCompSP27(para->getParD(level)->numberofthreads,
                               para->getParD(level)->omega,
                               para->getParD(level)->neighborX_SP,
                               para->getParD(level)->neighborY_SP,
                               para->getParD(level)->neighborZ_SP,
                               para->getParD(level)->d0SP.f[0],
                               para->getParD(level)->size_Mat_SP,
                               level,
                               para->getForcesDev(),
                               pm[i]->getPorosity(),
                               pm[i]->getDarcyLBM(),
                               pm[i]->getForchheimerLBM(),
                               pm[i]->getSizePM(),
                               pm[i]->getHostNodeIDsPM(),
                               para->getParD(level)->evenOrOdd);
        getLastCudaError("KernelPMCumOneCompSP27 execution failed");
    }
}

void collisionAdvectionDiffusion(Parameter* para, int level)
{
    if (para->getDiffMod() == 7)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // incompressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //KernelADincomp7(para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP,
        //                para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,
        //                para->getParD(level)->evenOrOdd);
        //getLastCudaError("KernelADincomp7 execution failed");

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // compressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //KernelThS7( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP,
        //            para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //            para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],      para->getParD(level)->size_Mat_SP,
        //            para->getParD(level)->evenOrOdd);
        //getLastCudaError("KernelThS7 execution failed");
    }
    else if (para->getDiffMod() == 27)
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // incompressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //KernelADincomp27( para->getParD(level)->numberofthreads,    para->getParD(level)->diffusivity,  para->getParD(level)->geoSP,
        //                    para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],     para->getParD(level)->size_Mat_SP,
        //                    para->getParD(level)->evenOrOdd);
        //getLastCudaError("KernelADincomp27 execution failed");

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // compressible
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //KernelThS27(para->getParD(level)->numberofthreads,
        //            para->getParD(level)->diffusivity,
        //            para->getParD(level)->geoSP,
        //            para->getParD(level)->neighborX_SP,
        //            para->getParD(level)->neighborY_SP,
        //            para->getParD(level)->neighborZ_SP,
        //            para->getParD(level)->d0SP.f[0],
        //            para->getParD(level)->d27.f[0],
        //            para->getParD(level)->size_Mat_SP,
        //            para->getParD(level)->evenOrOdd);
        //getLastCudaError("KernelThS27 execution failed");
    }
}

void prepareExchangeMultiGPU(Parameter *para, int level, int streamIndex)
{
    prepareExchangeCollDataXGPU27AllNodes(para, level, streamIndex);
    prepareExchangeCollDataYGPU27AllNodes(para, level, streamIndex);
    prepareExchangeCollDataZGPU27AllNodes(para, level, streamIndex);
}

void prepareExchangeMultiGPUAfterFtoC(Parameter *para, int level, int streamIndex)
{
    prepareExchangeCollDataXGPU27AfterFtoC(para, level, streamIndex);
    prepareExchangeCollDataYGPU27AfterFtoC(para, level, streamIndex);
    prepareExchangeCollDataZGPU27AfterFtoC(para, level, streamIndex);
}

void exchangeMultiGPU(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager, int level,
                      int streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AllNodes(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataYGPU27AllNodes(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataZGPU27AllNodes(para, comm, cudaManager, level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AllNodes(para, level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AllNodes(para, level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AllNodes(para, level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADYGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADZGPU27(para, comm, cudaManager, level);
    }

    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // 1D domain decomposition
    // exchangePostCollDataGPU27(para, comm, level);
}
void exchangeMultiGPU_noStreams_withPrepare(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager, int level, bool useReducedComm)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    if (useReducedComm) {
        // X
        prepareExchangeCollDataXGPU27AfterFtoC(para, level, -1);
        exchangeCollDataXGPU27AfterFtoC(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferXGPU27AfterFtoC(para, level, -1);
        // Y
        prepareExchangeCollDataYGPU27AfterFtoC(para, level, -1);
        exchangeCollDataYGPU27AfterFtoC(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferYGPU27AfterFtoC(para, level, -1);
        // Z
        prepareExchangeCollDataZGPU27AfterFtoC(para, level, -1);
        exchangeCollDataZGPU27AfterFtoC(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferZGPU27AfterFtoC(para, level, -1);
    } else {
        // X
        prepareExchangeCollDataXGPU27AllNodes(para, level, -1);
        exchangeCollDataXGPU27AllNodes(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferXGPU27AllNodes(para, level, -1);
        // Y
        prepareExchangeCollDataYGPU27AllNodes(para, level, -1);
        exchangeCollDataYGPU27AllNodes(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferYGPU27AllNodes(para, level, -1);
        // Z
        prepareExchangeCollDataZGPU27AllNodes(para, level, -1);
        exchangeCollDataZGPU27AllNodes(para, comm, cudaManager, level, -1);
        scatterNodesFromRecvBufferZGPU27AllNodes(para, level, -1);
    }

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADYGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADZGPU27(para, comm, cudaManager, level);
    }
}
void exchangeMultiGPUAfterFtoC(Parameter *para, vf::gpu::Communicator &comm, CudaMemoryManager *cudaManager, int level,
                               int streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataYGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataZGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);

    scatterNodesFromRecvBufferXGPU27AfterFtoC(para, level, streamIndex);
    scatterNodesFromRecvBufferYGPU27AfterFtoC(para, level, streamIndex);
    scatterNodesFromRecvBufferZGPU27AfterFtoC(para, level, streamIndex);

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition convection diffusion
    if (para->getDiffOn()) {
        if (para->getUseStreams())
            std::cout << "Warning: Cuda streams not yet implemented for convection diffusion" << std::endl;
        exchangePostCollDataADXGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADYGPU27(para, comm, cudaManager, level);
        exchangePostCollDataADZGPU27(para, comm, cudaManager, level);
    }
}

void UpdateGrid27::postCollisionBC(int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    // V E L O C I T Y (I N F L O W)
    this->cudaKernelManager->runVelocityBCKernel(level);

    //////////////////////////////////////////////////////////////////////////
    // N O - S L I P
    this->cudaKernelManager->runNoSlipBCKernel(level);

    //////////////////////////////////////////////////////////////////////////
    // S L I P
    this->cudaKernelManager->runSlipBCKernel(level);

    //////////////////////////////////////////////////////////////////////////
    // S T R E S S (wall model)
    this->cudaKernelManager->runStressWallModelKernel(level);

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    this->cudaKernelManager->runGeoBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    this->cudaKernelManager->runOutflowBCKernel(level);

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    this->cudaKernelManager->runPressureBCKernelPost(level);

    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    if (para->getDiffOn())
    {
        if (para->getDiffMod() == 7)
        {
            if (para->getParD(level)->QGeom.numberOfBCnodes > 0)
            {
                //QNoSlipADincompDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
                //                     para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,
                //                     para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0],
                //                     para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                //                     para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                //                     para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
                //getLastCudaError("QNoSlipADincompDev7 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADDev7( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
                         para->getParD(level)->d0SP.f[0],             para->getParD(level)->d7.f[0],      para->getParD(level)->Temp.temp,
                         para->getParD(level)->diffusivity,           para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0],
                         para->getParD(level)->Temp.kTemp,            para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                         para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                         para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
                getLastCudaError("QADDev27 execution failed");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if (para->getParD(level)->TempVel.kTemp > 0)
            {
                //QADVeloIncompDev7(  para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],            para->getParD(level)->TempVel.tempPulse,
                //                    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->TempVel.k,
                //                    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,    para->getParD(level)->TempVel.kTemp,
                //                    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                //                    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADVeloIncompDev7 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////

                //if (t<15580)//(t>500000 && t<515580)//(t>300000 && t<315580)
                //{
                //    QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                      para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],            para->getParD(level)->TempVel.tempPulse,
                //                      para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->TempVel.k,
                //                      para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,    para->getParD(level)->TempVel.kTemp,
                //                      para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,
                //                      para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev7 execution failed");
                //}
                //else
                //{
                //    QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                      para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],            para->getParD(level)->TempVel.temp,
                //                      para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->TempVel.k,
                //                      para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,    para->getParD(level)->TempVel.kTemp,
                //                      para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,
                //                      para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev7 execution failed");
                //}

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADVelDev7( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                            para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],            para->getParD(level)->TempVel.temp,
                            para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->TempVel.k,
                            para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,     para->getParD(level)->TempVel.kTemp,
                            para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                            para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVelDev27 execution failed");
            }

            if (para->getParD(level)->TempPress.kTemp > 0)
            {
                //QADPressIncompDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,                para->getParD(level)->ny,
                //                    para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],            para->getParD(level)->TempPress.temp,
                //                    para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,        para->getParD(level)->TempPress.k,
                //                    para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp, para->getParD(level)->TempPress.kTemp,
                //                    para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,
                //                    para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADPressIncompDev7 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADPressDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,                para->getParD(level)->ny,
                              para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],            para->getParD(level)->TempPress.temp,
                              para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,        para->getParD(level)->TempPress.k,
                              para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,   para->getParD(level)->TempPress.kTemp,
                              para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                              para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                getLastCudaError("QADPressDev27 execution failed");
            }
        }
        else if (para->getDiffMod() == 27)
        {
            if (para->getParD(level)->QGeom.numberOfBCnodes > 0)
            {
                //QNoSlipADincompDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
                //                     para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,
                //                     para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0],
                //                     para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                //                     para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                //                     para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
                //getLastCudaError("QNoSlipADincompDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADBBDev27(para->getParD(level)->numberofthreads,      para->getParD(level)->nx,           para->getParD(level)->ny,
                           para->getParD(level)->d0SP.f[0],            para->getParD(level)->d27.f[0],     para->getParD(level)->Temp.temp,
                           para->getParD(level)->diffusivity,          para->getParD(level)->Temp.k,       para->getParD(level)->QGeom.q27[0],
                           para->getParD(level)->Temp.kTemp,           para->getParD(level)->Temp.kTemp,   para->getParD(level)->omega,
                           para->getParD(level)->neighborX_SP,         para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP,          para->getParD(level)->evenOrOdd);
                getLastCudaError("QADBBDev27 execution failed");
            }

            if (para->getParD(level)->TempVel.kTemp > 0)
            {
                QADVeloIncompDev27( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                                    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],        para->getParD(level)->TempVel.temp,
                                    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->TempVel.k,
                                    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,
                                    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,
                                    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVeloIncompDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////

                //if (t>500000 && t<515580)//(t>300000 && t<315580)
                //{
                //    QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                 para->getParD(level)->ny,
                //                       para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],         para->getParD(level)->TempVel.tempPulse,
                //                       para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,     para->getParD(level)->TempVel.k,
                //                       para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,
                //                       para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP,
                //                       para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev27 execution failed");
                //}
                //else
                //{
                //    QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                       para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],        para->getParD(level)->TempVel.temp,
                //                       para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,    para->getParD(level)->TempVel.k,
                //                       para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,  para->getParD(level)->TempVel.kTemp,
                //                       para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,
                //                       para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,    para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev27 execution failed");
                //}

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                             para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],            para->getParD(level)->TempVel.tempPulse,
                             para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->Qinflow.k,
                             para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->numberOfVeloBCnodes,         para->getParD(level)->numberOfVeloBCnodes,
                             para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP,
                             para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVelDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////

                //if (t<1000)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
                //{
                //    QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                   para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],            para->getParD(level)->TempVel.tempPulse,
                //                   para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->Qinflow.k,
                //                   para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->numberOfVeloBCnodes,         para->getParD(level)->numberOfVeloBCnodes,
                //                   para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP,
                //                   para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVelDev27 execution failed");
                //}
                //else
                //{
                //    QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,                para->getParD(level)->ny,
                //                para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],            para->getParD(level)->TempVel.temp,
                //                para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,        para->getParD(level)->Qinflow.k,
                //                para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->numberOfVeloBCnodes,         para->getParD(level)->numberOfVeloBCnodes,
                //                para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                //                para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVelDev27 execution failed");
                //}
            }

            if (para->getParD(level)->TempPress.kTemp > 0)
            {
                //QADPressIncompDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,                para->getParD(level)->ny,
                //                     para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],            para->getParD(level)->TempPress.temp,
                //                     para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,        para->getParD(level)->TempPress.k,
                //                     para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,
                //                     para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                //                     para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADPressIncompDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADPressDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,                para->getParD(level)->ny,
                               para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],            para->getParD(level)->TempPress.temp,
                               para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,        para->getParD(level)->TempPress.k,
                               para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,
                               para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,
                               para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,        para->getParD(level)->evenOrOdd);
                getLastCudaError("QADPressDev27 execution failed");

            }
        }
    }
}

void swapBetweenEvenAndOddTimestep(Parameter* para, int level)
{
    if (para->getParD(level)->evenOrOdd==true)  para->getParD(level)->evenOrOdd=false;
    else                                        para->getParD(level)->evenOrOdd=true;
}

void calcMacroscopicQuantities(Parameter* para, int level)
{
    CalcMacCompSP27(para->getParD(level)->vx_SP,
                    para->getParD(level)->vy_SP,
                    para->getParD(level)->vz_SP,
                    para->getParD(level)->rho_SP,
                    para->getParD(level)->press_SP,
                    para->getParD(level)->geoSP,
                    para->getParD(level)->neighborX_SP,
                    para->getParD(level)->neighborY_SP,
                    para->getParD(level)->neighborZ_SP,
                    para->getParD(level)->size_Mat_SP,
                    para->getParD(level)->numberofthreads,
                    para->getParD(level)->d0SP.f[0],
                    para->getParD(level)->evenOrOdd);
    getLastCudaError("CalcMacSP27 execution failed");
}

void preCollisionBC(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    // I N F L O W
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->numberOfVeloBCnodes > 0)
    {
        // TODO: https://git.rz.tu-bs.de/irmb/VirtualFluids_dev/-/issues/29
        //if (  myid == 0)
        //{
        //    VelSchlaffer27(para->getParD(level)->numberofthreads, t,
        //                   para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.Vz,
        //                   para->getParD(level)->Qinflow.deltaVz, para->getParD(level)->Qinflow.k,
        //                   para->getParD(level)->Qinflow.kN,      para->getParD(level)->numberOfVeloBCnodes,
        //                   para->getParD(level)->omega,           para->getParD(level)->neighborX_SP,
        //                   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
        //                   para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //    getLastCudaError("VelSchlaffer27 execution failed");
        //}
        //////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        //QVelDevIncompHighNu27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                      para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //                      para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0],
        //                      para->getParD(level)->numberOfVeloBCnodes,        para->getParD(level)->numberOfVeloBCnodes,     para->getParD(level)->omega,
        //                      para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                      para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDevComp27 execution failed");
        //////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        //QVelDevCompHighNu27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                    para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //                    para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0],
        //                    para->getParD(level)->numberOfVeloBCnodes,        para->getParD(level)->numberOfVeloBCnodes,     para->getParD(level)->omega,
        //                    para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                    para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDevComp27 execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->QGeom.numberOfBCnodes > 0)
    {
        if (para->getCalcDragLift())
        {
            //Drag and Lift Part II
            DragLiftPreD27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->QGeom.k,
                           para->getParD(level)->QGeom.q27[0],
                           para->getParD(level)->QGeom.numberOfBCnodes,
                           para->getParD(level)->DragPreX,
                           para->getParD(level)->DragPreY,
                           para->getParD(level)->DragPreZ,
                           para->getParD(level)->neighborX_SP,
                           para->getParD(level)->neighborY_SP,
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP,
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
            getLastCudaError("DragLift27 execution failed");
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of Drag and Lift
            ////////////////////////////////////////////////////////////////////////////////
            calcDragLift(para, cudaManager, level);
            ////////////////////////////////////////////////////////////////////////////////
        }

        if (para->getCalcCp())
        {
            ////////////////////////////////////////////////////////////////////////////////
            //Calculation of cp
            ////////////////////////////////////////////////////////////////////////////////

            if(t > para->getTStartOut())
            {
                ////////////////////////////////////////////////////////////////////////////////
                CalcCPtop27(para->getParD(level)->d0SP.f[0],
                            para->getParD(level)->cpTopIndex,
                            para->getParD(level)->numberOfPointsCpTop,
                            para->getParD(level)->cpPressTop,
                            para->getParD(level)->neighborX_SP,
                            para->getParD(level)->neighborY_SP,
                            para->getParD(level)->neighborZ_SP,
                            para->getParD(level)->size_Mat_SP,
                            para->getParD(level)->evenOrOdd,
                            para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(para->getParD(level)->d0SP.f[0],
                               para->getParD(level)->cpBottomIndex,
                               para->getParD(level)->numberOfPointsCpBottom,
                               para->getParD(level)->cpPressBottom,
                               para->getParD(level)->neighborX_SP,
                               para->getParD(level)->neighborY_SP,
                               para->getParD(level)->neighborZ_SP,
                               para->getParD(level)->size_Mat_SP,
                               para->getParD(level)->evenOrOdd,
                               para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                CalcCPbottom27(para->getParD(level)->d0SP.f[0],
                               para->getParD(level)->cpBottom2Index,
                               para->getParD(level)->numberOfPointsCpBottom2,
                               para->getParD(level)->cpPressBottom2,
                               para->getParD(level)->neighborX_SP,
                               para->getParD(level)->neighborY_SP,
                               para->getParD(level)->neighborZ_SP,
                               para->getParD(level)->size_Mat_SP,
                               para->getParD(level)->evenOrOdd,
                               para->getParD(level)->numberofthreads);
                //////////////////////////////////////////////////////////////////////////////////
                calcCp(para, cudaManager, level);
            }
        }


        ////////////////////////////////////////////////////////////////////////////////
        // high viscosity incompressible
        //QDevIncompHighNu27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //                    para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0],
        //                    para->getParD(level)->QGeom.numberOfBCnodes,              para->getParD(level)->QGeom.numberOfBCnodes,     para->getParD(level)->omega,
        //                    para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                    para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDevComp27 (Geom) execution failed");
        //////////////////////////////////////////////////////////////////////////////////
        // high viscosity compressible
        //QDevCompHighNu27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //                  para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0],
        //                  para->getParD(level)->QGeom.numberOfBCnodes,              para->getParD(level)->QGeom.numberOfBCnodes,     para->getParD(level)->omega,
        //                  para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                  para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDevComp27 (Geom) execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->QPress.numberOfBCnodes > 0)
    {
        QPressNoRhoDev27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
                        para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
                        para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.numberOfBCnodes,     para->getParD(level)->omega,
                        para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,  para->getParD(level)->neighborZ_SP,
                        para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        getLastCudaError("QPressNoRhoDev27 execution failed");

        //QPressDevEQZ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
        //               para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
        //               para->getParD(level)->QPress.kN,       para->getParD(level)->kDistTestRE.f[0],
        //               para->getParD(level)->QPress.numberOfBCnodes,       para->getParD(level)->omega,
        //               para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(0)->neighborZ_SP,
        //               para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressDevEQZ27 execution failed");

        //QInflowScaleByPressDev27( para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
        //                          para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
        //                          para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.numberOfBCnodes,    para->getParD(0)->omega,
        //                          para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(0)->neighborZ_SP,
        //                          para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QInflowScaleByPressDev27 execution failed");

        ////////////////////////////////////////////////////////////////////////////////
        //press NEQ incompressible
        //QPressDevIncompNEQ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
        //                     para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
        //                     para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.numberOfBCnodes,    para->getParD(level)->omega,
        //                     para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                     para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressDevIncompNEQ27 execution failed");
        //////////////////////////////////////////////////////////////////////////////////
        //press NEQ compressible
        // QPressDevNEQ27( para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
        //                 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
        //                 para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.numberOfBCnodes,    para->getParD(level)->omega,
        //                 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        // getLastCudaError("QPressDevNEQ27 execution failed");

    }

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->numberOfOutflowBCnodes > 0)
    {
        //QPressNoRhoDev27( para->getParD(level)->numberofthreads, para->getParD(level)->Qoutflow.RhoBC,
        //                  para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qoutflow.k,
        //                  para->getParD(level)->Qoutflow.kN,     para->getParD(level)->Qoutflow.numberOfBCnodes,    para->getParD(level)->omega,
        //                  para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
        //                  para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressNoRhoDev27 execution failed");

        //if (  myid == numprocs - 1)
        //PressSchlaffer27( para->getParD(level)->numberofthreads,  para->getParD(level)->Qoutflow.RhoBC,
        //                  para->getParD(level)->d0SP.f[0],        para->getParD(level)->Qoutflow.Vx,
        //                  para->getParD(level)->Qoutflow.Vy,      para->getParD(level)->Qoutflow.Vz,
        //                  para->getParD(level)->Qoutflow.deltaVz, para->getParD(level)->Qoutflow.k,
        //                  para->getParD(level)->Qoutflow.kN,      para->getParD(level)->numberOfOutflowBCnodes,
        //                  para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,
        //                  para->getParD(level)->neighborY_SP,     para->getParD(level)->neighborZ_SP,
        //                  para->getParD(level)->size_Mat_SP,      para->getParD(level)->evenOrOdd);
        //getLastCudaError("PressSchlaffer27 execution failed");
    }

    //////////////////////////////////////////////////////////////////////////////////
    ////only for a round off error test
    //para->cudaCopyTestREtoHost(0,para->getParH(0)->QPress.numberOfBCnodes);
    //printRE(para, t);
    //////////////////////////////////////////////////////////////////////////////////


}

void fineToCoarse(Parameter* para, int level)
{
    //ScaleFC_comp_D3Q27F3(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],     para->getParD(level)->g6.g[0],
    //                     para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                     para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                     para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                     para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                     para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                     para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                     para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                     para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC_comp_D3Q27F3 execution failed");

    //ScaleFC_0817_comp_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                     para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                     para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                     para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                     para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                     para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                     para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                     para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                     para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC_0817_comp_27 execution failed");

    //ScaleFC_RhoSq_3rdMom_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                                para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                                para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                                para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                                para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                                para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                                para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                                para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                                para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC_RhoSq_3rdMom_comp_27 execution failed");

    ScaleFC_RhoSq_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                            para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
                            para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
                            para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
                            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
                            para->getParD(level)->offFC,          CU_STREAM_LEGACY);
    getLastCudaError("ScaleFC27_RhoSq_comp execution failed");

    //ScaleFC_AA2016_comp_27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                        para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                        para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                        para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                        para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                        para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                        para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                        para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                        para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC_AA2016_comp_27 execution failed");



    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //ScaleFC27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
    //            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP, para->getParD(level+1)->neighborZ_SP,
    //            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
    //            para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //            para->getParD(level)->K_FC,           para->getParD(level)->omega,          para->getParD(level+1)->omega,
    //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny,
    //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
    //getLastCudaError("ScaleFC27 execution failed");

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

    // ScaleFC_Fix_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                      para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                      para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                      para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                      para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                      para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                      para->getParD(level)->offFC);
    // getLastCudaError("ScaleFC27 execution failed");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // incompressible
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

    //ScaleFC_NSPress_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                   para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                   para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                   para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                   para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
    //                   para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                   para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                   para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                   para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC27 execution failed");


    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
    {
        if (para->getDiffMod() == 7)
        {
            //ScaleFCThS7(   para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
            //               para->getParD(level)->d7.f[0],        para->getParD(level+1)->d7.f[0],
            //               para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
            //               para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
            //               para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
            //               para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF,
            //               para->getParD(level)->K_FC,
            //               para->getParD(level)->vis,            para->getParD(level)->diffusivity,     para->getParD(level)->numberofthreads);
            //getLastCudaError("ScaleFCTh7 execution failed");

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
            ScaleFCThS27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
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

}

void fineToCoarseWithStream(Parameter *para, int level, uint *iCellFCC, uint *iCellFCF, uint k_FC, int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager()->getStream(streamIndex);

    ScaleFC_RhoSq_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                            iCellFCC,                             iCellFCF,
                            k_FC,                                 para->getParD(level)->omega,           para->getParD(level + 1)->omega,
                            para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
                            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
                            para->getParD(level)->offFC,          stream);
    getLastCudaError("ScaleFC27_RhoSq_comp execution failed");

    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn()) {
        printf("fineToCoarseWithStream Advection Diffusion not implemented"); // TODO
    }
}

void coarseToFine(Parameter* para, int level)
{
    //ScaleCF_comp_D3Q27F3(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],     para->getParD(level+1)->g6.g[0],
    //                     para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                     para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                     para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                     para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                     para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                     para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                     para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                     para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF_comp_D3Q27F3 execution failed");

    //ScaleCF_0817_comp_27(para->getParD(level)->d0SP.f[0],      para->getParD(level + 1)->d0SP.f[0],
    //                     para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                     para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                     para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                     para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                     para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                     para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                     para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                     para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF_0817_comp_27 execution failed");

    //ScaleCF_RhoSq_3rdMom_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                              para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                              para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                              para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                              para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                              para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                              para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                              para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                              para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF_RhoSq_3rdMom_comp_27 execution failed");

    ScaleCF_RhoSq_comp_27(para->getParD(level)->d0SP.f[0],        para->getParD(level + 1)->d0SP.f[0],
                          para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP,     para->getParD(level)->neighborZ_SP,
                          para->getParD(level + 1)->neighborX_SP, para->getParD(level + 1)->neighborY_SP, para->getParD(level + 1)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP,      para->getParD(level + 1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
                          para->getParD(level)->intCF.ICellCFC,   para->getParD(level)->intCF.ICellCFF,
                          para->getParD(level)->K_CF,             para->getParD(level)->omega,            para->getParD(level + 1)->omega,
                          para->getParD(level)->vis,              para->getParD(level)->nx,               para->getParD(level)->ny,
                          para->getParD(level + 1)->nx,           para->getParD(level + 1)->ny,           para->getParD(level)->numberofthreads,
                          para->getParD(level)->offCF,            cudaStreamLegacy);
    getLastCudaError("ScaleCF27_RhoSq_comp execution failed");

    //ScaleCF_AA2016_comp_27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                        para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                        para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                        para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                        para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                        para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                        para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                        para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                        para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF_AA2016_comp_27 execution failed");



    //////////////////////////////////////////////////////////////////////////
    // D E P R E C A T E D
    //////////////////////////////////////////////////////////////////////////

    //ScaleCF27(  para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
    //            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP, para->getParD(level+1)->neighborZ_SP,
    //            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
    //            para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //            para->getParD(level)->K_CF,           para->getParD(level)->omega,          para->getParD(level+1)->omega,
    //            para->getParD(level)->vis,            para->getParD(level)->nx,             para->getParD(level)->ny,
    //            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,           para->getParD(level)->gridNX);
    //getLastCudaError("ScaleCF27 execution failed");

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

    // ScaleCF_Fix_comp_27(    para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                      para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                      para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                      para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                      para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                      para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                      para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                      para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                      para->getParD(level)->offCF);
    // getLastCudaError("ScaleCF27 execution failed");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // incompressible
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

    //ScaleCF_NSPress_27(para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
    //                   para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
    //                   para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
    //                   para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
    //                   para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
    //                   para->getParD(level)->K_CF,           para->getParD(level)->omega,           para->getParD(level+1)->omega,
    //                   para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny,
    //                   para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
    //                   para->getParD(level)->offCF);
    //getLastCudaError("ScaleCF27 execution failed");


    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    //////////////////////////////////////////////////////////////////////////

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
        }
        else if (para->getDiffMod() == 27)
        {
            ScaleCFThS27( para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],
                          para->getParD(level)->d27.f[0],       para->getParD(level+1)->d27.f[0],
                          para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP,
                          para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
                          para->getParD(level)->intCF.ICellCFC, para->getParD(level)->intCF.ICellCFF,
                          para->getParD(level)->K_CF,
                          para->getParD(level)->vis,            para->getParD(level+1)->diffusivity,   para->getParD(level)->numberofthreads,
                          para->getParD(level)->offCF);
            getLastCudaError("ScaleCFTh27 execution failed");
        }
    }

}

void coarseToFineWithStream(Parameter *para, int level, uint *iCellCFC, uint *iCellCFF, uint k_CF, OffCF &offCF,
                            int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager()->getStream(streamIndex);

    ScaleCF_RhoSq_comp_27(para->getParD(level)->d0SP.f[0],        para->getParD(level + 1)->d0SP.f[0],
                          para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP,     para->getParD(level)->neighborZ_SP,
                          para->getParD(level + 1)->neighborX_SP, para->getParD(level + 1)->neighborY_SP, para->getParD(level + 1)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP,      para->getParD(level + 1)->size_Mat_SP,  para->getParD(level)->evenOrOdd,
                          iCellCFC,                               iCellCFF,
                          k_CF,                                   para->getParD(level)->omega,            para->getParD(level + 1)->omega,
                          para->getParD(level)->vis,              para->getParD(level)->nx,               para->getParD(level)->ny,
                          para->getParD(level + 1)->nx,           para->getParD(level + 1)->ny,           para->getParD(level)->numberofthreads,
                          offCF,                                  stream);
    getLastCudaError("ScaleCF27_RhoSq_comp execution failed");

    if (para->getDiffOn()) {
        printf("CoarseToFineWithStream Advection Diffusion not implemented"); // TODO
    }
}


UpdateGrid27::UpdateGrid27(SPtr<Parameter> para, vf::gpu::Communicator &comm, SPtr<CudaMemoryManager> cudaManager,
                           std::vector<std::shared_ptr<PorousMedia>> &pm, std::vector<SPtr<Kernel>> &kernels)
    : para(para), comm(comm), cudaManager(cudaManager), pm(pm), kernels(kernels)
{
    chooseFunctionForCollisionAndExchange();
    chooseFunctionForRefinementAndExchange();
    this->cudaKernelManager = CudaKernelManager::make(para);
}


void UpdateGrid27::chooseFunctionForCollisionAndExchange()
{
    std::cout << "Function used for collisionAndExchange: ";
    if (para->getUseStreams() && para->getNumprocs() > 1 && para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = &UpdateGrid27::collisionAndExchange_streams;
        std::cout << "collisionAndExchange_streams()" << std::endl;

    } else if (para->getUseStreams() && !para->getKernelNeedsFluidNodeIndicesToRun()) {
        std::cout << "Cuda Streams can only be used with kernels which run using fluidNodesIndices." << std::endl;

    } else if (para->getUseStreams() && para->getNumprocs() <= 1) {
        std::cout << "Cuda Streams can only be used with multiple MPI processes." << std::endl;

    } else if (!para->getUseStreams() && para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = &UpdateGrid27::collisionAndExchange_noStreams_indexKernel;
        std::cout << "collisionAndExchange_noStreams_indexKernel()" << std::endl;

    } else if (!para->getUseStreams() && !para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = &UpdateGrid27::collisionAndExchange_noStreams_oldKernel;
        std::cout << "collisionAndExchange_noStreams_oldKernel()" << std::endl;

    } else {
        std::cout << "Invalid Configuration for collision and exchange" << std::endl;
    }
}

void UpdateGrid27::chooseFunctionForRefinementAndExchange()
{
    std::cout << "Function used for refinementAndExchange: ";
    if (para->getMaxLevel() == 0) {
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_noRefinementAndExchange;
        std::cout << "only one level - no function needed." << std::endl;

    } else if (para->getNumprocs() == 1) {
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_noExchange;
        std::cout << "refinementAndExchange_noExchange()" << std::endl;

    } else if (para->getNumprocs() > 1 && para->getUseStreams() && para->useReducedCommunicationAfterFtoC) {
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_streams_onlyExchangeInterface;
        std::cout << "refinementAndExchange_streams_onlyExchangeInterface()" << std::endl;

    } else if(para->getNumprocs() > 1 && para->getUseStreams() && !para->useReducedCommunicationAfterFtoC){
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_streams_completeExchange;
        std::cout << "refinementAndExchange_streams_completeExchange()" << std::endl;

    } else if (para->getNumprocs() > 1 && !para->getUseStreams() && para->useReducedCommunicationAfterFtoC) {
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_noStreams_onlyExchangeInterface;
        std::cout << "refinementAndExchange_noStreams_onlyExchangeInterface()" << std::endl;

    } else {
        this->refinementAndExchange = &UpdateGrid27::refinementAndExchange_noStreams_completeExchange;
        std::cout << "refinementAndExchange_noStreams_completeExchange()" << std::endl;
    }
}

void interactWithActuators(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> actuator: para->getActuators() )
    {
        actuator->interact(para, cudaManager, level, t);
    }
}

void interactWithProbes(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)
{
    for( SPtr<PreCollisionInteractor> probe: para->getProbes() )
    {
        probe->interact(para, cudaManager, level, t);
    }
}

void calcTurbulentViscosity(Parameter* para, int level)
{
    if(para->getUseAMD())
        calcTurbulentViscosityAMD(para, level);
}