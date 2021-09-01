#include "Communication/ExchangeData27.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataXGPU27(Parameter *para, int level, int streamIndex) 
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
        GetSendFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->sendProcessNeighborX[i].f[0],
                           para->getParD(level)->sendProcessNeighborX[i].index,
                           para->getParD(level)->sendProcessNeighborX[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads,
                           stream);    
}

void exchangeCollDataXGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager, int level,
                                int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
        cudaManager->cudaCopyProcessNeighborXFsDH(level, i, streamIndex);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborX[i].f[0],
                            para->getParH(level)->recvProcessNeighborX[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborX[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborX[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborX[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborX[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait for memcopy device to host to finish before sending data
    if (para->getUseStreams())
        cudaStreamSynchronize(stream); 
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborX[i].f[0],
                          para->getParH(level)->sendProcessNeighborX[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborX[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborXFsHD(level, i, streamIndex);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->recvProcessNeighborX[i].f[0],
                           para->getParD(level)->recvProcessNeighborX[i].index,
                           para->getParD(level)->recvProcessNeighborX[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads,
                           stream);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataYGPU27(Parameter *para, int level, int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);   
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
        GetSendFsPostDev27(para->getParD(level)->d0SP.f[0], 
                           para->getParD(level)->sendProcessNeighborY[i].f[0],
                           para->getParD(level)->sendProcessNeighborY[i].index,
                           para->getParD(level)->sendProcessNeighborY[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP,
                           para->getParD(level)->neighborZ_SP, 
                           para->getParD(level)->size_Mat_SP,
                           para->getParD(level)->evenOrOdd, 
                           para->getParD(level)->numberofthreads, 
                           stream);
}

void exchangeCollDataYGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager, int level,
                                int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
        cudaManager->cudaCopyProcessNeighborYFsDH(level, i, streamIndex);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborY[i].f[0],
                            para->getParH(level)->recvProcessNeighborY[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborY[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborY[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait for memcopy device to host to finish before sending data
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborY[i].f[0],
                          para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborYFsHD(level, i, streamIndex);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->recvProcessNeighborY[i].f[0],
                           para->getParD(level)->recvProcessNeighborY[i].index,
                           para->getParD(level)->recvProcessNeighborY[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads,
                           stream);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void prepareExchangeCollDataZGPU27(Parameter *para, int level, int streamIndex) {
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);   
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
        GetSendFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->sendProcessNeighborZ[i].f[0],
                           para->getParD(level)->sendProcessNeighborZ[i].index,
                           para->getParD(level)->sendProcessNeighborZ[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads,
                           stream);
} 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataZGPU27(Parameter *para, vf::gpu::Communicator *comm, CudaMemoryManager *cudaManager, int level,
                                int streamIndex)
{
    cudaStream_t stream = (streamIndex == -1) ? CU_STREAM_LEGACY : para->getStreamManager().getStream(streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
        cudaManager->cudaCopyProcessNeighborZFsDH(level, i, streamIndex);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborZ[i].f[0],
                            para->getParH(level)->recvProcessNeighborZ[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborZ[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborZ[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // wait for memcopy device to host to finish before sending data
    if (para->getUseStreams())
        cudaStreamSynchronize(stream);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborZ[i].f[0],
                          para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborZFsHD(level, i, streamIndex);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->recvProcessNeighborZ[i].f[0],
                           para->getParD(level)->recvProcessNeighborZ[i].index,
                           para->getParD(level)->recvProcessNeighborZ[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads,
                           stream);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



















////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//1D domain decomposition
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighbors(level, "send")); i++)
    {
        //////////////////////////////////////////////////////////////////////////
        GetSendFsPreDev27(para->getParD(level)->d0SP.f[0],
                          para->getParD(level)->sendProcessNeighbor[i].f[0],
                          para->getParD(level)->sendProcessNeighbor[i].index,
                          para->getParD(level)->sendProcessNeighbor[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborFsDH(level, i);
        //////////////////////////////////////////////////////////////////////////
        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighbor[i].f[0], 
                            para->getParH(level)->sendProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighbor[i].f[0],
                            para->getParH(level)->recvProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->sendProcessNeighbor[i].rankNeighbor);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->d0SP.f[0],
                          para->getParD(level)->recvProcessNeighbor[i].f[0],
                          para->getParD(level)->recvProcessNeighbor[i].index,
                          para->getParD(level)->recvProcessNeighbor[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighbors(level, "send")); i++)
    {
        //////////////////////////////////////////////////////////////////////////
        GetSendFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->sendProcessNeighbor[i].f[0],
                           para->getParD(level)->sendProcessNeighbor[i].index,
                           para->getParD(level)->sendProcessNeighbor[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborFsDH(level, i);
        //////////////////////////////////////////////////////////////////////////
        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighbor[i].f[0], 
                            para->getParH(level)->sendProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighbor[i].f[0],
                            para->getParH(level)->recvProcessNeighbor[i].numberOfFs,
                            para->getParH(level)->sendProcessNeighbor[i].rankNeighbor);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d0SP.f[0],
                           para->getParD(level)->recvProcessNeighbor[i].f[0],
                           para->getParD(level)->recvProcessNeighbor[i].index,
                           para->getParD(level)->recvProcessNeighbor[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////3D domain decomposition
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// X
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePreCollDataXGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->sendProcessNeighborX[i].f[0],
//                          para->getParD(level)->sendProcessNeighborX[i].index,
//                          para->getParD(level)->sendProcessNeighborX[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborX[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborX[i].f[0],
//                            para->getParH(level)->recvProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborX[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->recvProcessNeighborX[i].f[0],
//                          para->getParD(level)->recvProcessNeighborX[i].index,
//                          para->getParD(level)->recvProcessNeighborX[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataXGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->sendProcessNeighborX[i].f[0],
//                            para->getParD(level)->sendProcessNeighborX[i].index,
//                            para->getParD(level)->sendProcessNeighborX[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborX[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborX[i].f[0],
//                            para->getParH(level)->recvProcessNeighborX[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborX[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborXFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->recvProcessNeighborX[i].f[0],
//                            para->getParD(level)->recvProcessNeighborX[i].index,
//                            para->getParD(level)->recvProcessNeighborX[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePreCollDataYGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->sendProcessNeighborY[i].f[0],
//                          para->getParD(level)->sendProcessNeighborY[i].index,
//                          para->getParD(level)->sendProcessNeighborY[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborY[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborY[i].f[0],
//                            para->getParH(level)->recvProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->recvProcessNeighborY[i].f[0],
//                          para->getParD(level)->recvProcessNeighborY[i].index,
//                          para->getParD(level)->recvProcessNeighborY[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataYGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->sendProcessNeighborY[i].f[0],
//                            para->getParD(level)->sendProcessNeighborY[i].index,
//                            para->getParD(level)->sendProcessNeighborY[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborY[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborY[i].f[0],
//                            para->getParH(level)->recvProcessNeighborY[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborY[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborYFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->recvProcessNeighborY[i].f[0],
//                            para->getParD(level)->recvProcessNeighborY[i].index,
//                            para->getParD(level)->recvProcessNeighborY[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePreCollDataZGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->sendProcessNeighborZ[i].f[0],
//                          para->getParD(level)->sendProcessNeighborZ[i].index,
//                          para->getParD(level)->sendProcessNeighborZ[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborZ[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParH(level)->recvProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPreDev27(para->getParD(level)->d0SP.f[0],
//                          para->getParD(level)->recvProcessNeighborZ[i].f[0],
//                          para->getParD(level)->recvProcessNeighborZ[i].index,
//                          para->getParD(level)->recvProcessNeighborZ[i].numberOfNodes,
//                          para->getParD(level)->neighborX_SP, 
//                          para->getParD(level)->neighborY_SP, 
//                          para->getParD(level)->neighborZ_SP,
//                          para->getParD(level)->size_Mat_SP, 
//                          para->getParD(level)->evenOrOdd,
//                          para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void exchangePostCollDataZGPU27(Parameter* para, vf::gpu::Communicator* comm, int level)
//{
//    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
//    {
//        //////////////////////////////////////////////////////////////////////////
//        GetSendFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->sendProcessNeighborZ[i].f[0],
//                            para->getParD(level)->sendProcessNeighborZ[i].index,
//                            para->getParD(level)->sendProcessNeighborZ[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsDH(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        comm->exchngDataGPU(para->getParH(level)->sendProcessNeighborZ[i].f[0], 
//                            para->getParH(level)->sendProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParH(level)->recvProcessNeighborZ[i].numberOfFs,
//                            para->getParH(level)->sendProcessNeighborZ[i].rankNeighbor);
//        //////////////////////////////////////////////////////////////////////////
//        para->cudaCopyProcessNeighborZFsHD(level, i);
//        //////////////////////////////////////////////////////////////////////////
//        SetRecvFsPostDev27( para->getParD(level)->d0SP.f[0],
//                            para->getParD(level)->recvProcessNeighborZ[i].f[0],
//                            para->getParD(level)->recvProcessNeighborZ[i].index,
//                            para->getParD(level)->recvProcessNeighborZ[i].numberOfNodes,
//                            para->getParD(level)->neighborX_SP, 
//                            para->getParD(level)->neighborY_SP, 
//                            para->getParD(level)->neighborZ_SP,
//                            para->getParD(level)->size_Mat_SP, 
//                            para->getParD(level)->evenOrOdd,
//                            para->getParD(level)->numberofthreads);
//        //////////////////////////////////////////////////////////////////////////
//    }
//}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


















































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3D domain decomposition convection diffusion
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADXGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->sendProcessNeighborADX[i].f[0],
                          para->getParD(level)->sendProcessNeighborADX[i].index,
                          para->getParD(level)->sendProcessNeighborADX[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADXFsDH(level, i);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADX[i].f[0],
                            para->getParH(level)->recvProcessNeighborADX[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADX[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
                          para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADXFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->recvProcessNeighborADX[i].f[0],
                          para->getParD(level)->recvProcessNeighborADX[i].index,
                          para->getParD(level)->recvProcessNeighborADX[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADXGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->sendProcessNeighborADX[i].f[0],
                           para->getParD(level)->sendProcessNeighborADX[i].index,
                           para->getParD(level)->sendProcessNeighborADX[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADXFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADX[i].f[0],
                            para->getParH(level)->recvProcessNeighborADX[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADX[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADX[i].f[0],
                          para->getParH(level)->sendProcessNeighborADX[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADX[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADXFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->recvProcessNeighborADX[i].f[0],
                           para->getParD(level)->recvProcessNeighborADX[i].index,
                           para->getParD(level)->recvProcessNeighborADX[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADYGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->sendProcessNeighborADY[i].f[0],
                          para->getParD(level)->sendProcessNeighborADY[i].index,
                          para->getParD(level)->sendProcessNeighborADY[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADYFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADY[i].f[0],
                            para->getParH(level)->recvProcessNeighborADY[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADY[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
                          para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADYFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->recvProcessNeighborADY[i].f[0],
                          para->getParD(level)->recvProcessNeighborADY[i].index,
                          para->getParD(level)->recvProcessNeighborADY[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADYGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->sendProcessNeighborADY[i].f[0],
                           para->getParD(level)->sendProcessNeighborADY[i].index,
                           para->getParD(level)->sendProcessNeighborADY[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADYFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADY[i].f[0],
                            para->getParH(level)->recvProcessNeighborADY[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADY[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADY[i].f[0],
                          para->getParH(level)->sendProcessNeighborADY[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADY[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADYFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->recvProcessNeighborADY[i].f[0],
                           para->getParD(level)->recvProcessNeighborADY[i].index,
                           para->getParD(level)->recvProcessNeighborADY[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePreCollDataADZGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        GetSendFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParD(level)->sendProcessNeighborADZ[i].index,
                          para->getParD(level)->sendProcessNeighborADZ[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADZ[i].f[0],
                            para->getParH(level)->recvProcessNeighborADZ[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADZ[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPreDev27(para->getParD(level)->d27.f[0],
                          para->getParD(level)->recvProcessNeighborADZ[i].f[0],
                          para->getParD(level)->recvProcessNeighborADZ[i].index,
                          para->getParD(level)->recvProcessNeighborADZ[i].numberOfNodes,
                          para->getParD(level)->neighborX_SP, 
                          para->getParD(level)->neighborY_SP, 
                          para->getParD(level)->neighborZ_SP,
                          para->getParD(level)->size_Mat_SP, 
                          para->getParD(level)->evenOrOdd,
                          para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangePostCollDataADZGPU27(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        GetSendFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->sendProcessNeighborADZ[i].f[0],
                           para->getParD(level)->sendProcessNeighborADZ[i].index,
                           para->getParD(level)->sendProcessNeighborADZ[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborADZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->nbRecvDataGPU(para->getParH(level)->recvProcessNeighborADZ[i].f[0],
                            para->getParH(level)->recvProcessNeighborADZ[i].numberOfFs,
                            para->getParH(level)->recvProcessNeighborADZ[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////start non blocking MPI send
    //for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    //{
    //    comm->nbSendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
    //                        para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
    //                        para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Waitall
    //if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    //{
    //    comm->waitallGPU();
    //}
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->sendDataGPU(para->getParH(level)->sendProcessNeighborADZ[i].f[0],
                          para->getParH(level)->sendProcessNeighborADZ[i].numberOfFs,
                          para->getParH(level)->sendProcessNeighborADZ[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborADZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        SetRecvFsPostDev27(para->getParD(level)->d27.f[0],
                           para->getParD(level)->recvProcessNeighborADZ[i].f[0],
                           para->getParD(level)->recvProcessNeighborADZ[i].index,
                           para->getParD(level)->recvProcessNeighborADZ[i].numberOfNodes,
                           para->getParD(level)->neighborX_SP, 
                           para->getParD(level)->neighborY_SP, 
                           para->getParD(level)->neighborZ_SP,
                           para->getParD(level)->size_Mat_SP, 
                           para->getParD(level)->evenOrOdd,
                           para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

















































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3D domain decomposition F3 - K18/K20
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// X
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3XGPU(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3X[i].g[0],
            para->getParD(level)->sendProcessNeighborF3X[i].index,
            para->getParD(level)->sendProcessNeighborF3X[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborF3XFsDH(level, i);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->nbRecvDataGPU(
            para->getParH(level)->recvProcessNeighborF3X[i].g[0],
            para->getParH(level)->recvProcessNeighborF3X[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3X[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->sendDataGPU(
            para->getParH(level)->sendProcessNeighborF3X[i].g[0],
            para->getParH(level)->sendProcessNeighborF3X[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3X[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsX(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborF3XFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3X[i].g[0],
            para->getParD(level)->recvProcessNeighborF3X[i].index,
            para->getParD(level)->recvProcessNeighborF3X[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Y
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3YGPU(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3Y[i].g[0],
            para->getParD(level)->sendProcessNeighborF3Y[i].index,
            para->getParD(level)->sendProcessNeighborF3Y[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborF3YFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->nbRecvDataGPU(
            para->getParH(level)->recvProcessNeighborF3Y[i].g[0],
            para->getParH(level)->recvProcessNeighborF3Y[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3Y[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->sendDataGPU(
            para->getParH(level)->sendProcessNeighborF3Y[i].g[0],
            para->getParH(level)->sendProcessNeighborF3Y[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3Y[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsY(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborF3YFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3Y[i].g[0],
            para->getParD(level)->recvProcessNeighborF3Y[i].index,
            para->getParD(level)->recvProcessNeighborF3Y[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Z
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void exchangeCollDataF3ZGPU(Parameter* para, vf::gpu::Communicator* comm, CudaMemoryManager* cudaManager, int level)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Device to Host
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        getSendGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->sendProcessNeighborF3Z[i].g[0],
            para->getParD(level)->sendProcessNeighborF3Z[i].index,
            para->getParD(level)->sendProcessNeighborF3Z[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
        //////////////////////////////////////////////////////////////////////////
        cudaManager->cudaCopyProcessNeighborF3ZFsDH(level, i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start non blocking MPI receive
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->nbRecvDataGPU(
            para->getParH(level)->recvProcessNeighborF3Z[i].g[0],
            para->getParH(level)->recvProcessNeighborF3Z[i].numberOfGs,
            para->getParH(level)->recvProcessNeighborF3Z[i].rankNeighbor);
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //start blocking MPI send
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->sendDataGPU(
            para->getParH(level)->sendProcessNeighborF3Z[i].g[0],
            para->getParH(level)->sendProcessNeighborF3Z[i].numberOfGs,
            para->getParH(level)->sendProcessNeighborF3Z[i].rankNeighbor);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Wait
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        comm->waitGPU(i);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //reset the request array
    if (0 < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")))
    {
        comm->resetRequest();
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //copy Host to Device
    for (unsigned int i = 0; i < (unsigned int)(para->getNumberOfProcessNeighborsZ(level, "send")); i++)
    {
        cudaManager->cudaCopyProcessNeighborF3ZFsHD(level, i);
        //////////////////////////////////////////////////////////////////////////
        setRecvGsDevF3(
            para->getParD(level)->g6.g[0],
            para->getParD(level)->recvProcessNeighborF3Z[i].g[0],
            para->getParD(level)->recvProcessNeighborF3Z[i].index,
            para->getParD(level)->recvProcessNeighborF3Z[i].numberOfNodes,
            para->getParD(level)->neighborX_SP,
            para->getParD(level)->neighborY_SP,
            para->getParD(level)->neighborZ_SP,
            para->getParD(level)->size_Mat_SP,
            para->getParD(level)->evenOrOdd,
            para->getParD(level)->numberofthreads);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




















