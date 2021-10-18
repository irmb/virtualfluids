#include "Calculation/UpdateGrid27.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "Calculation/DragLift.h"
#include "Calculation/Cp.h"
//#include "Utilities/StringUtil.hpp"
//#include "Output/UnstructuredGridWriter.hpp"
#include "Communication/ExchangeData27.h"
#include "Kernel/Kernel.h"
#include "Parameter/CudaStreamManager.h"

void UpdateGrid27::updateGrid27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager,
                  std::vector<std::shared_ptr<PorousMedia>> &pm, int level, unsigned int t,
                  std::vector<SPtr<Kernel>> &kernels)
{
    //////////////////////////////////////////////////////////////////////////

    if (level != para->getFine()) {
        updateGrid27(para, comm, cudaManager, pm, level + 1, t, kernels);
        updateGrid27(para, comm, cudaManager, pm, level + 1, t, kernels);
    }

    //////////////////////////////////////////////////////////////////////////

    this->collisionAndExchange(para, pm, level, t, kernels, comm, cudaManager);

    //////////////////////////////////////////////////////////////////////////

    postCollisionBC(para, level, t);

    //////////////////////////////////////////////////////////////////////////

    swapBetweenEvenAndOddTimestep(para, level);

	//////////////////////////////////////////////////////////////////////////

	if (para->getUseWale())
		calcMacroscopicQuantities(para, level);

	//////////////////////////////////////////////////////////////////////////

    preCollisionBC(para, cudaManager, level, t);

    //////////////////////////////////////////////////////////////////////////
    if( level != para->getFine() )
    {
        this->refinementAndExchange(para, level, comm, cudaManager);
    }
}

void refinementAndExchange_streams(Parameter *para, int level, vf::gpu::Communicator *comm,
                                   CudaMemoryManager *cudaManager)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex   = para->getStreamManager()->getBulkStreamIndex();

    fineToCoarseWithStream(para, level, para->getParD(level)->intFCBorder.ICellFCC,
                           para->getParD(level)->intFCBorder.ICellFCF, para->getParD(level)->intFCBorder.kFC,
                           borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    prepareExchangeMultiGPUAfterFtoC(para, level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernel
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    fineToCoarseWithStream(para, level, para->getParD(level)->intFCBulk.ICellFCC,
                           para->getParD(level)->intFCBulk.ICellFCF, para->getParD(level)->intFCBulk.kFC,
                           bulkStreamIndex);

    exchangeMultiGPUAfterFtoC(para, comm, cudaManager, level, borderStreamIndex);
    coarseToFine(para, level);
}

void refinementAndExchange_noStreams_onlyExchangeInterface(Parameter *para, int level, vf::gpu::Communicator *comm,
                                                           CudaMemoryManager *cudaManager)
{
    fineToCoarse(para, level);

    prepareExchangeMultiGPUAfterFtoC(para, level, -1);
    exchangeMultiGPUAfterFtoC(para, comm, cudaManager, level, -1);

    coarseToFine(para, level);
}

void refinementAndExchange_noStreams_completeExchange(Parameter *para, int level, vf::gpu::Communicator *comm,
                                                      CudaMemoryManager *cudaManager)
{
    fineToCoarse(para, level);

    prepareExchangeMultiGPU(para, level, -1);
    exchangeMultiGPU(para, comm, cudaManager, level, -1);

    coarseToFine(para, level);
}

void refinementAndExchange_noExchange(Parameter *para, int level, vf::gpu::Communicator *comm,
                                      CudaMemoryManager *cudaManager)
{
    fineToCoarse(para, level);
    coarseToFine(para, level);
}

void collisionAndExchange_noStreams_indexKernel(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm,
                                                       int level, unsigned int t, std::vector<SPtr<Kernel>> &kernels,
                                                       vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager)
{
    collisionUsingIndex(para, pm, level, t, kernels, para->getParD(level)->fluidNodeIndices,
                            para->getParD(level)->numberOfFluidNodes, -1);
    prepareExchangeMultiGPU(para, level, -1);
    exchangeMultiGPU(para, comm, cudaManager, level, -1);
}

void collisionAndExchange_noStreams_oldKernel(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm,
                                                 int level, unsigned int t, std::vector<SPtr<Kernel>> &kernels,
                                                 vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager)
{
    collision(para, pm, level, t, kernels);
    prepareExchangeMultiGPU(para, level, -1);
    exchangeMultiGPU(para, comm, cudaManager, level, -1);
}

void collisionAndExchange_streams(Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level,
                                     unsigned int t, std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager)
{
    int borderStreamIndex = para->getStreamManager()->getBorderStreamIndex();
    int bulkStreamIndex   = para->getStreamManager()->getBulkStreamIndex();
    // launch border kernel
    collisionUsingIndex(para, pm, level, t, kernels, para->getParD(level)->fluidNodeIndicesBorder,
                        para->getParD(level)->numberOffluidNodesBorder, borderStreamIndex);

    // prepare exchange and trigger bulk kernel when finished
    prepareExchangeMultiGPU(para, level, borderStreamIndex);
    if (para->getUseStreams())
        para->getStreamManager()->triggerStartBulkKernel(borderStreamIndex);

    // launch bulk kernel
    para->getStreamManager()->waitOnStartBulkKernelEvent(bulkStreamIndex);
    collisionUsingIndex(para, pm, level, t, kernels, para->getParD(level)->fluidNodeIndices,
                        para->getParD(level)->numberOfFluidNodes, bulkStreamIndex);

    exchangeMultiGPU(para, comm, cudaManager, level, borderStreamIndex);
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
        //		            para->getParD(level)->neighborX_SP,       para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //		            para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],     para->getParD(level)->size_Mat_SP,  
        //		            para->getParD(level)->evenOrOdd); 
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

void exchangeMultiGPU(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager, int level,
                      int streamIndex)
{
    // St. Lenz: exchange for post-collision data and pre-collision data are identical!

    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AllNodes(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataYGPU27AllNodes(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataZGPU27AllNodes(para, comm, cudaManager, level, streamIndex);

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
void exchangeMultiGPUAfterFtoC(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager, int level,
                               int streamIndex)
{
    //////////////////////////////////////////////////////////////////////////
    // 3D domain decomposition
    exchangeCollDataXGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataYGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);
    exchangeCollDataZGPU27AfterFtoC(para, comm, cudaManager, level, streamIndex);

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

void postCollisionBC(Parameter* para, int level, unsigned int t)
{
    //////////////////////////////////////////////////////////////////////////
    // I N F L O W
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->kInflowQ > 0)
    {
        //QVelDev27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //           para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //           para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
        //           para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
        //           para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //           para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDev27 execution failed");
        
        QVelDevComp27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
                       para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
                       para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
                       para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
                       para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                       para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        getLastCudaError("QVelDevComp27 execution failed");

        //QVelDevCompZeroPress27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,             para->getParD(level)->ny,
        //                       para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,     para->getParD(level)->Qinflow.Vz,
        //                       para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,      para->getParD(level)->Qinflow.q27[0],
        //                       para->getParD(level)->kInflowQ,        para->getParD(level)->Qinflow.kArray, para->getParD(level)->omega,
        //                       para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
        //                       para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDevCompZeroPress27 execution failed");

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////

        //QVelDevice1h27( para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //                para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
        //                para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,          
        //                para->getPhi(),                        para->getAngularVelocity(),
        //                para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                para->getParD(level)->coordX_SP,       para->getParD(level)->coordY_SP,    para->getParD(level)->coordZ_SP,
        //                para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDev27 execution failed");

    }

    //////////////////////////////////////////////////////////////////////////
    // N O - S L I P
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->kQ > 0)
    {
        //QDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //	      para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
        //	      para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
        //	      para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //	      para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDev27 execution failed");
        
        //BBDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //         para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
        //         para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
        //         para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //         para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("BBDev27 (Wall) execution failed");

        //QDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //        para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
        //        para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
        //        para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //        para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDev27 (Wall) execution failed");

        QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
                   para->getParD(level)->d0SP.f[0],             para->getParD(level)->QWall.k,      para->getParD(level)->QWall.q27[0], 
                   para->getParD(level)->kQ,                    para->getParD(level)->kQ,           para->getParD(level)->omega,
                   para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                   para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        getLastCudaError("QDevComp27 (Wall) execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // S L I P
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->kSlipQ > 0)
    {
        //QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
        //            para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
        //            para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //            para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("Slip27 execution failed");

        QSlipDevComp27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QSlip.k,
                        para->getParD(level)->QSlip.q27[0],    para->getParD(level)->kSlipQ,       para->getParD(level)->omega,
                        para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                        para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        getLastCudaError("QSlipDev27 execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // G E O M E T R Y
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->QGeom.kQ > 0)
    {
		if (para->getCalcDragLift())
		{
			//Drag and Lift Part I
			DragLiftPostD27(para->getParD(level)->d0SP.f[0], 
			                para->getParD(level)->QGeom.k, 
			                para->getParD(level)->QGeom.q27[0],
			                para->getParD(level)->QGeom.kQ, 
			                para->getParD(level)->DragPostX,
			                para->getParD(level)->DragPostY,
			                para->getParD(level)->DragPostZ,
			                para->getParD(level)->neighborX_SP,
			                para->getParD(level)->neighborY_SP,
			                para->getParD(level)->neighborZ_SP,
			                para->getParD(level)->size_Mat_SP, 
			                para->getParD(level)->evenOrOdd,
			                para->getParD(level)->numberofthreads);
			getLastCudaError("DragLift27 execution failed"); 
		}

        //BBDev27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //         para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
        //         para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //         para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //         para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("BBDev27 (Wall) execution failed");

        //QDev27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //		 para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
        //		 para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //		 para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //		 para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDev27 (Geom) execution failed");

        //QVelDev27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //          para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
        //          para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
        //          para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //          para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //          para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDev27 execution failed");

        //QDevComp27(para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //           para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
        //           para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //           para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //           para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDevComp27 (Geom) execution failed");

        QVelDevComp27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
                      para->getParD(level)->QGeom.Vx,        para->getParD(level)->QGeom.Vy,     para->getParD(level)->QGeom.Vz,
                      para->getParD(level)->d0SP.f[0],       para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
                      para->getParD(level)->QGeom.kQ,        para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
                      para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                      para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        getLastCudaError("QVelDevComp27 execution failed");

        //QVelDevCompZeroPress27(	para->getParD(0)->numberofthreads, para->getParD(0)->nx,           para->getParD(0)->ny,
		//						para->getParD(0)->QGeom.Vx,        para->getParD(0)->QGeom.Vy,     para->getParD(0)->QGeom.Vz,
		//						para->getParD(0)->d0SP.f[0],       para->getParD(0)->QGeom.k,      para->getParD(0)->QGeom.q27[0], 
		//						para->getParD(0)->QGeom.kQ,        para->getParD(0)->QGeom.kQ,     para->getParD(0)->omega,
		//						para->getParD(0)->neighborX_SP,    para->getParD(0)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//						para->getParD(0)->size_Mat_SP,     para->getParD(0)->evenOrOdd);
		//getLastCudaError("QVelDevCompZeroPress27 execution failed");

        //QDev3rdMomentsComp27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
        //                      para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
        //                      para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //                      para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                      para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
        //getLastCudaError("QDev3rdMomentsComp27 (Geom) execution failed");

        //QSlipDev27( para->getParD(level)->numberofthreads, para->getParD(level)->d0SP.f[0],    para->getParD(level)->QGeom.k,
        //            para->getParD(level)->QGeom.q27[0],    para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
        //            para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //            para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("Slip27 execution failed");

        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////
        // the GridGenerator does currently not provide normals!
        
        //QSlipGeomDevComp27(para->getParD(level)->numberofthreads,     para->getParD(level)->d0SP.f[0],           para->getParD(level)->QGeom.k,
        //                   para->getParD(level)->QGeom.q27[0],        para->getParD(level)->QGeom.kQ,            para->getParD(level)->omega,
        //                   para->getParD(level)->QGeomNormalX.q27[0], para->getParD(level)->QGeomNormalY.q27[0], para->getParD(level)->QGeomNormalZ.q27[0],
        //                   para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,        para->getParD(level)->neighborZ_SP,
        //                   para->getParD(level)->size_Mat_SP,         para->getParD(level)->evenOrOdd);
        //getLastCudaError("QSlipGeomDev27 execution failed");

        //QSlipNormDevComp27(para->getParD(level)->numberofthreads,     para->getParD(level)->d0SP.f[0],           para->getParD(level)->QGeom.k,
        //                   para->getParD(level)->QGeom.q27[0],        para->getParD(level)->QGeom.kQ,            para->getParD(level)->omega,
        //                   para->getParD(level)->QGeomNormalX.q27[0], para->getParD(level)->QGeomNormalY.q27[0], para->getParD(level)->QGeomNormalZ.q27[0],
        //                   para->getParD(level)->neighborX_SP,        para->getParD(level)->neighborY_SP,        para->getParD(level)->neighborZ_SP,
        //                   para->getParD(level)->size_Mat_SP,         para->getParD(level)->evenOrOdd);
        //getLastCudaError("QSlipGeomDev27 execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->kOutflowQ > 0)
    {
        //////////////////////////////////////////////////////////////////////////
        // D E P R E C A T E D
        //////////////////////////////////////////////////////////////////////////

        //QPressDevFixBackflow27( para->getParD(level)->numberofthreads,       RhoBCOutflowD,
        //                        para->getParD(level)->d0SP.f[0],              QoutflowD.k, kOutflowQ,             para->getParD(level)->omega,
        //                        para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                        para->getParD(level)->size_Mat_SP,  para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressDev27 execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    //////////////////////////////////////////////////////////////////////////

    if (para->getParD(level)->kPressQ > 0)
    {
        //QPressDev27_IntBB(  para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
        //					para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,       para->getParD(level)->QPress.q27[0], 
        //					para->getParD(level)->QPress.kQ,       para->getParD(level)->QPress.kQ,      para->getParD(level)->omega,
        //					para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
        //					para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressDev27_IntBB fine execution failed");
    }

    //////////////////////////////////////////////////////////////////////////
    // A D V E C T I O N    D I F F U S I O N
    //////////////////////////////////////////////////////////////////////////

    if (para->getDiffOn())
    {
        if (para->getDiffMod() == 7)
        {
            if (para->getParD(level)->QGeom.kQ > 0)
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
                //QADVeloIncompDev7(  para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                //                    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
                //                    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
                //                    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
                //                    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	    para->getParD(level)->neighborY_SP, 
                //                    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADVeloIncompDev7 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////

                //if (t<15580)//(t>500000 && t<515580)//(t>300000 && t<315580)
                //{
                //    QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                //                      para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.tempPulse, 
                //                      para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
                //                      para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
                //                      para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
                //                      para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev7 execution failed");
                //}
                //else
                //{
                //    QADVeloIncompDev7(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                //                      para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
                //                      para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
                //                      para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,	para->getParD(level)->TempVel.kTemp,  
                //                      para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
                //                      para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev7 execution failed");
                //}

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////
                
                QADVelDev7( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                            para->getParD(level)->d0SP.f[0],          para->getParD(level)->d7.f[0],			para->getParD(level)->TempVel.temp, 
                            para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
                            para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,     para->getParD(level)->TempVel.kTemp,  
                            para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
                            para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVelDev27 execution failed");
            }

            if (para->getParD(level)->TempPress.kTemp > 0)
            {
                //QADPressIncompDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
                //                    para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
                //                    para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
                //                    para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp, para->getParD(level)->TempPress.kTemp,  
                //                    para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
                //                    para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADPressIncompDev7 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////

                QADPressDev7( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
                              para->getParD(level)->d0SP.f[0],        para->getParD(level)->d7.f[0],			para->getParD(level)->TempPress.temp, 
                              para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
                              para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,   para->getParD(level)->TempPress.kTemp,  
                              para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
                              para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                getLastCudaError("QADPressDev27 execution failed");
            }
        } 
        else if (para->getDiffMod() == 27)
        {
            if (para->getParD(level)->QGeom.kQ > 0)
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
                QADVeloIncompDev27( para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
            	                    para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
            	                    para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->TempVel.k,
            	                    para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
            	                    para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
            	                    para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVeloIncompDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////
                
                //if (t>500000 && t<515580)//(t>300000 && t<315580)
                //{
                //    QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				 para->getParD(level)->ny,
                //                       para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		 para->getParD(level)->TempVel.tempPulse, 
                //                       para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	 para->getParD(level)->TempVel.k,
                //                       para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,   para->getParD(level)->TempVel.kTemp,  
                //                       para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	 para->getParD(level)->neighborY_SP, 
                //                       para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	 para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev27 execution failed");
                //}
                //else
                //{
                //    QADVeloIncompDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                //                       para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		para->getParD(level)->TempVel.temp, 
                //                       para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,	para->getParD(level)->TempVel.k,
                //                       para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->TempVel.kTemp,  para->getParD(level)->TempVel.kTemp,  
                //                       para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	para->getParD(level)->neighborY_SP, 
                //                       para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,	para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVeloIncompDev27 execution failed");
                //}

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////
                
                QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                 	        para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		    para->getParD(level)->TempVel.tempPulse, 
                 	        para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
                 	        para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,         para->getParD(level)->kInflowQ,  
                 	        para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP, 
                 	        para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                getLastCudaError("QADVelDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // D E P R E C A T E D
                //////////////////////////////////////////////////////////////////////////

                //if (t<1000)//(t>100000 && t<103895)//(t>1600000 && t<1662317)//(t>500000 && t<515580)//(t<1000)//(t<15580)//(t>400000 && t<415580)//
                //{
                //    QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                // 	              para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		    para->getParD(level)->TempVel.tempPulse, 
                // 	              para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
                // 	              para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,         para->getParD(level)->kInflowQ,  
                // 	              para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,     para->getParD(level)->neighborY_SP, 
                // 	              para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVelDev27 execution failed");
                //}
                //else
                //{
                //    QADVelDev27(para->getParD(level)->numberofthreads,    para->getParD(level)->nx,				para->getParD(level)->ny,
                //                para->getParD(level)->d0SP.f[0],          para->getParD(level)->d27.f[0],		    para->getParD(level)->TempVel.temp, 
                //                para->getParD(level)->TempVel.velo,       para->getParD(level)->diffusivity,		para->getParD(level)->Qinflow.k,
                //                para->getParD(level)->Qinflow.q27[0],     para->getParD(level)->kInflowQ,         para->getParD(level)->kInflowQ,  
                //                para->getParD(level)->omega,              para->getParD(level)->neighborX_SP,	    para->getParD(level)->neighborY_SP, 
                //                para->getParD(level)->neighborZ_SP,       para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //    getLastCudaError("QADVelDev27 execution failed");
                //}
            }

            if (para->getParD(level)->TempPress.kTemp > 0)
            {
                //QADPressIncompDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
                //                     para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
                //                     para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
                //                     para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
                //                     para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
                //                     para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
                //getLastCudaError("QADPressIncompDev27 execution failed");

                //////////////////////////////////////////////////////////////////////////
                // C O M P R E S S I B L E
                //////////////////////////////////////////////////////////////////////////
                
                QADPressDev27( para->getParD(level)->numberofthreads,  para->getParD(level)->nx,				para->getParD(level)->ny,
                               para->getParD(level)->d0SP.f[0],        para->getParD(level)->d27.f[0],			para->getParD(level)->TempPress.temp, 
                               para->getParD(level)->TempPress.velo,   para->getParD(level)->diffusivity,		para->getParD(level)->TempPress.k,
                               para->getParD(level)->QPress.q27[0],    para->getParD(level)->TempPress.kTemp,  para->getParD(level)->TempPress.kTemp,  
                               para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,		para->getParD(level)->neighborY_SP, 
                               para->getParD(level)->neighborZ_SP,     para->getParD(level)->size_Mat_SP,		para->getParD(level)->evenOrOdd);
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

	if (para->getParD(level)->kInflowQ > 0)
	{
		//if (  myid == 0)
		//{
		//    VelSchlaffer27(para->getParD(level)->numberofthreads, t,
		//                   para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.Vz, 
		//                   para->getParD(level)->Qinflow.deltaVz, para->getParD(level)->Qinflow.k,  
		//                   para->getParD(level)->Qinflow.kN,      para->getParD(level)->kInflowQ, 
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
		//                      para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
		//                      para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		//                      para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QVelDevComp27 execution failed");
		//////////////////////////////////////////////////////////////////////////////
		// high viscosity compressible
		//QVelDevCompHighNu27(para->getParD(level)->numberofthreads, para->getParD(level)->nx,           para->getParD(level)->ny,
        //                    para->getParD(level)->Qinflow.Vx,      para->getParD(level)->Qinflow.Vy,   para->getParD(level)->Qinflow.Vz,
        //                    para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qinflow.k,    para->getParD(level)->Qinflow.q27[0], 
        //                    para->getParD(level)->kInflowQ,        para->getParD(level)->kInflowQ,     para->getParD(level)->omega,
        //                    para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                    para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QVelDevComp27 execution failed");
	}

	//////////////////////////////////////////////////////////////////////////
	// G E O M E T R Y
	//////////////////////////////////////////////////////////////////////////

	if (para->getParD(level)->QGeom.kQ > 0)
	{
		if (para->getCalcDragLift())
		{
			//Drag and Lift Part II
			DragLiftPreD27(para->getParD(level)->d0SP.f[0], 
			               para->getParD(level)->QGeom.k, 
			               para->getParD(level)->QGeom.q27[0],
			               para->getParD(level)->QGeom.kQ, 
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
		//                    para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		//                    para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		//                    para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		//getLastCudaError("QDevComp27 (Geom) execution failed");
		//////////////////////////////////////////////////////////////////////////////////
		// high viscosity compressible
		//QDevCompHighNu27( para->getParD(level)->numberofthreads,       para->getParD(level)->nx,           para->getParD(level)->ny,
		//                  para->getParD(level)->d0SP.f[0],             para->getParD(level)->QGeom.k,      para->getParD(level)->QGeom.q27[0], 
		//                  para->getParD(level)->QGeom.kQ,              para->getParD(level)->QGeom.kQ,     para->getParD(level)->omega,
		//                  para->getParD(level)->neighborX_SP,          para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
		//                  para->getParD(level)->size_Mat_SP,           para->getParD(level)->evenOrOdd);
		//getLastCudaError("QDevComp27 (Geom) execution failed");
	}

    //////////////////////////////////////////////////////////////////////////
    // P R E S S U R E
    //////////////////////////////////////////////////////////////////////////

	if (para->getParD(level)->QPress.kQ > 0)
	{
		//QPressNoRhoDev27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC,
		//                 para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,
		//                 para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,     para->getParD(level)->omega,
		//                 para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,  para->getParD(level)->neighborZ_SP,
		//                 para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QPressNoRhoDev27 execution failed");

		//QPressDevEQZ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
		//               para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
		//               para->getParD(level)->QPress.kN,       para->getParD(level)->kDistTestRE.f[0],       
		//               para->getParD(level)->QPress.kQ,       para->getParD(level)->omega,
		//               para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//               para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QPressDevEQZ27 execution failed");

		//QInflowScaleByPressDev27( para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
		//                          para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
		//                          para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,    para->getParD(0)->omega,
		//                          para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(0)->neighborZ_SP,
		//                          para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QInflowScaleByPressDev27 execution failed");

		////////////////////////////////////////////////////////////////////////////////
        //press NEQ incompressible
        //QPressDevIncompNEQ27(para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
        //                     para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
        //                     para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,    para->getParD(level)->omega,
        //                     para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
        //                     para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        //getLastCudaError("QPressDevIncompNEQ27 execution failed");
        //////////////////////////////////////////////////////////////////////////////////
        //press NEQ compressible
        QPressDevNEQ27( para->getParD(level)->numberofthreads, para->getParD(level)->QPress.RhoBC, 
                        para->getParD(level)->d0SP.f[0],       para->getParD(level)->QPress.k,  
                        para->getParD(level)->QPress.kN,       para->getParD(level)->QPress.kQ,    para->getParD(level)->omega,
                        para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP, para->getParD(level)->neighborZ_SP,
                        para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
        getLastCudaError("QPressDevNEQ27 execution failed");

	}

	//////////////////////////////////////////////////////////////////////////
    // O U T F L O W
    //////////////////////////////////////////////////////////////////////////

	if (para->getParD(level)->kOutflowQ > 0)
	{
		//QPressNoRhoDev27( para->getParD(level)->numberofthreads, para->getParD(level)->Qoutflow.RhoBC, 
		//                  para->getParD(level)->d0SP.f[0],       para->getParD(level)->Qoutflow.k,  
		//                  para->getParD(level)->Qoutflow.kN,     para->getParD(level)->Qoutflow.kQ,    para->getParD(level)->omega,
		//                  para->getParD(level)->neighborX_SP,    para->getParD(level)->neighborY_SP,   para->getParD(level)->neighborZ_SP,
		//                  para->getParD(level)->size_Mat_SP,     para->getParD(level)->evenOrOdd);
		//getLastCudaError("QPressNoRhoDev27 execution failed");

		//if (  myid == numprocs - 1)  
		//PressSchlaffer27( para->getParD(level)->numberofthreads,  para->getParD(level)->Qoutflow.RhoBC,
		//                  para->getParD(level)->d0SP.f[0],        para->getParD(level)->Qoutflow.Vx, 
		//                  para->getParD(level)->Qoutflow.Vy,      para->getParD(level)->Qoutflow.Vz, 
		//                  para->getParD(level)->Qoutflow.deltaVz, para->getParD(level)->Qoutflow.k,  
		//                  para->getParD(level)->Qoutflow.kN,      para->getParD(level)->kOutflowQ,                      
		//                  para->getParD(level)->omega,            para->getParD(level)->neighborX_SP,    
		//                  para->getParD(level)->neighborY_SP,     para->getParD(level)->neighborZ_SP,
		//                  para->getParD(level)->size_Mat_SP,      para->getParD(level)->evenOrOdd);
		//getLastCudaError("PressSchlaffer27 execution failed");
	}

    //////////////////////////////////////////////////////////////////////////////////
    ////only for a round off error test
    //para->cudaCopyTestREtoHost(0,para->getParH(0)->QPress.kQ);
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

	//ScaleFC_RhoSq_3rdMom_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
	//	                            para->getParD(level)->neighborX_SP,   para->getParD(level)->neighborY_SP,    para->getParD(level)->neighborZ_SP, 
	//	                            para->getParD(level+1)->neighborX_SP, para->getParD(level+1)->neighborY_SP,  para->getParD(level+1)->neighborZ_SP, 
	//	                            para->getParD(level)->size_Mat_SP,    para->getParD(level+1)->size_Mat_SP,   para->getParD(level)->evenOrOdd,
	//	                            para->getParD(level)->intFC.ICellFCC, para->getParD(level)->intFC.ICellFCF, 
	//	                            para->getParD(level)->K_FC,           para->getParD(level)->omega,           para->getParD(level+1)->omega, 
	//	                            para->getParD(level)->vis,            para->getParD(level)->nx,              para->getParD(level)->ny, 
	//	                            para->getParD(level+1)->nx,           para->getParD(level+1)->ny,            para->getParD(level)->numberofthreads,
	//	                            para->getParD(level)->offFC);
    //getLastCudaError("ScaleFC_RhoSq_3rdMom_comp_27 execution failed");

	ScaleFC_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
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

	// ScaleFC_Fix_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
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

    ScaleFC_RhoSq_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0], 
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
        if (para->getDiffMod() == 7) {
            // TODO
            printf("fineToCoarseWithStream Advection Diffusion not implemented");
            //ScaleFCThSMG7(para->getParD(level)->d0SP.f[0], para->getParD(level + 1)->d0SP.f[0],
            //              para->getParD(level)->d7.f[0], para->getParD(level + 1)->d7.f[0],
            //              para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP,
            //              para->getParD(level)->neighborZ_SP, para->getParD(level + 1)->neighborX_SP,
            //              para->getParD(level + 1)->neighborY_SP, para->getParD(level + 1)->neighborZ_SP,
            //              para->getParD(level)->size_Mat_SP, para->getParD(level + 1)->size_Mat_SP,
            //              para->getParD(level)->evenOrOdd, para->getParD(level)->intFC.ICellFCC,
            //              para->getParD(level)->intFC.ICellFCF, para->getParD(level)->K_FC, para->getParD(level)->vis,
            //              para->getParD(level)->diffusivity, para->getParD(level)->numberofthreads,
            //              para->getParD(level)->offFC);
            //getLastCudaError("ScaleFCTh7 execution failed");
        } else if (para->getDiffMod() == 27) {
            // TODO
            printf("fineToCoarseWithStream Advection Diffusion not implemented");
            //ScaleFCThS27(para->getParD(level)->d0SP.f[0], para->getParD(level + 1)->d0SP.f[0],
            //             para->getParD(level)->d27.f[0], para->getParD(level + 1)->d27.f[0],
            //             para->getParD(level)->neighborX_SP, para->getParD(level)->neighborY_SP,
            //             para->getParD(level)->neighborZ_SP, para->getParD(level + 1)->neighborX_SP,
            //             para->getParD(level + 1)->neighborY_SP, para->getParD(level + 1)->neighborZ_SP,
            //             para->getParD(level)->size_Mat_SP, para->getParD(level + 1)->size_Mat_SP,
            //             para->getParD(level)->evenOrOdd, para->getParD(level)->intFC.ICellFCC,
            //             para->getParD(level)->intFC.ICellFCF, para->getParD(level)->K_FC, para->getParD(level)->vis,
            //             para->getParD(level)->diffusivity, para->getParD(level)->numberofthreads,
            //             para->getParD(level)->offFC);
            //getLastCudaError("ScaleFCTh27 execution failed");
        }
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

	//ScaleCF_RhoSq_3rdMom_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
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
                          para->getParD(level)->offCF);
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
    
    // ScaleCF_Fix_comp_27(	para->getParD(level)->d0SP.f[0],      para->getParD(level+1)->d0SP.f[0],                
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

UpdateGrid27::UpdateGrid27() = default;
UpdateGrid27::~UpdateGrid27(){}

UpdateGrid27::UpdateGrid27(Parameter *para) { 
    chooseFunctionForCollisionAndExchange(para);
    chooseFunctionForRefinementAndExchange(para);
}


void UpdateGrid27::chooseFunctionForCollisionAndExchange(Parameter *para)
{
    std::cout << "Function used for collisionAndExchange: ";
    if (para->getUseStreams() && para->getNumprocs() > 1 && para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = [](Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level,
                                        unsigned int t, std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm,
                                        CudaMemoryManager *cudaManager) {
            collisionAndExchange_streams(para, pm, level, t, kernels, comm, cudaManager);
        };
        std::cout << "collisionAndExchange_streams()" << std::endl;

    } else if (para->getUseStreams() && !para->getKernelNeedsFluidNodeIndicesToRun()) {
        std::cout << "Cuda Streams can only be used with kernels which run using fluidNodesIndices." << std::endl;

    } else if (para->getUseStreams() && para->getNumprocs() <= 1) {
        std::cout << "Cuda Streams can only be with multiple MPI processes." << std::endl;
    
    } else if (!para->getUseStreams() && para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = [](Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level,
                                        unsigned int t, std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm,
                                        CudaMemoryManager *cudaManager) {
            collisionAndExchange_noStreams_indexKernel(para, pm, level, t, kernels, comm, cudaManager);
        };
        std::cout << "collisionAndExchange_noStreams_indexKernel()" << std::endl;
    
    } else if (!para->getUseStreams() && !para->getKernelNeedsFluidNodeIndicesToRun()) {
        this->collisionAndExchange = [](Parameter *para, std::vector<std::shared_ptr<PorousMedia>> &pm, int level,
                                        unsigned int t, std::vector<SPtr<Kernel>> &kernels, vf::gpu::Communicator *comm,
                                        CudaMemoryManager *cudaManager) {
            collisionAndExchange_noStreams_oldKernel(para, pm, level, t, kernels, comm, cudaManager);
        };
        std::cout << "collisionAndExchange_noStreams_oldKernel()" << std::endl;
    
    } else {
        std::cout << "Invalid Configuration for collision and exchange" << std::endl;
    }
}

void UpdateGrid27::chooseFunctionForRefinementAndExchange(Parameter *para)
{
    std::cout << "Function used for refinementAndExchange: ";
    if (para->getMaxLevel() == 0) {
        this->refinementAndExchange = [](Parameter *para, int level, vf::gpu::Communicator *comm,
                                         CudaMemoryManager *cudaManager) {};
        std::cout << "only one level - no function needed." << std::endl;
    } else if (para->getNumprocs() == 1){
        this->refinementAndExchange = [](Parameter *para, int level, vf::gpu::Communicator *comm,
                                         CudaMemoryManager *cudaManager) {
            refinementAndExchange_noExchange(para, level, comm, cudaManager);
        };
        std::cout << "refinementAndExchange_noExchange()" << std::endl;
    } else if (para->getUseStreams() && para->getNumprocs() > 1 && para->useReducedCommunicationAfterFtoC) {
        this->refinementAndExchange = [](Parameter *para, int level, vf::gpu::Communicator *comm,
                                         CudaMemoryManager *cudaManager) {
            refinementAndExchange_streams(para, level, comm, cudaManager);
        };
        std::cout << "refinementAndExchange_streams()" << std::endl;
    } else if (para->getNumprocs() > 1 && para->useReducedCommunicationAfterFtoC) {
        this->refinementAndExchange = [](Parameter *para, int level, vf::gpu::Communicator *comm,
                                         CudaMemoryManager *cudaManager) {
            refinementAndExchange_noStreams_onlyExchangeInterface(para, level, comm, cudaManager);
        };
        std::cout << "refinementAndExchange_noStreams_onlyExchangeInterface()" << std::endl;
    } else {
        this->refinementAndExchange = [](Parameter *para, int level, vf::gpu::Communicator *comm,
                                         CudaMemoryManager *cudaManager) {
            refinementAndExchange_noStreams_completeExchange(para, level, comm, cudaManager);
        };
        std::cout << "refinementAndExchange_noStreams_completeExchange()" << std::endl;
    }
}