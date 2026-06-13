//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
#include "F16CompressibleAdvectionDiffusion.h"

#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include <cuda_helper/CudaGrid.h>
#include <gpu/core/Utilities/RunAdvectionDiffusionCollision.cuh>
#include <lbm/advectionDiffusion/collision/F16AdvectionDiffusion.h>

using namespace vf::lbm::advection_diffusion;

namespace vf::gpu {

template <TurbulenceModel turbulenceModel>
std::shared_ptr<F16CompressibleAdvectionDiffusion<turbulenceModel>> F16CompressibleAdvectionDiffusion<turbulenceModel>::getNewInstance(std::shared_ptr<Parameter> para, int level)
{
    return std::shared_ptr<F16CompressibleAdvectionDiffusion>(new F16CompressibleAdvectionDiffusion(para, level));
}

template <TurbulenceModel turbulenceModel>
void F16CompressibleAdvectionDiffusion<turbulenceModel>::run()
{
    runOnIndices(para->getParD(level)->taggedFluidNodeIndices[CollisionTemplate::AllFeatures], para->getParD(level)->numberOfTaggedFluidNodes[CollisionTemplate::AllFeatures], CollisionTemplate::AllFeatures);
}

template <TurbulenceModel turbulenceModel>
void F16CompressibleAdvectionDiffusion<turbulenceModel>::runOnIndices(const uint* indices, uint size_indices,  CollisionTemplate /**/, CudaStreamIndex streamIdx)
{
    const auto parameter = vf::gpu::ad::getCollisionParameter(para->getParD(level).get(), para->getTurbulentPrandtlNumber(), indices, size_indices);

    auto collision = [] __device__(ADCollisionParameter & parameters) { return runF16AdvectionDiffusion(parameters); };

    cudaStream_t stream = para->getStreamManager()->getStream(streamIdx);
    const vf::cuda::CudaGrid grid(para->getParD(level)->numberofthreads, parameter.numberOfFluidNodes);

    vf::gpu::ad::runCollisionAdvectionDiffusion<decltype(collision), turbulenceModel><<<grid.grid, grid.threads, 0, stream>>>(collision, parameter);
    getLastCudaError("F16CompressibleAdvectionDiffusion execution failed");
}

template <TurbulenceModel turbulenceModel>
F16CompressibleAdvectionDiffusion<turbulenceModel>::F16CompressibleAdvectionDiffusion(std::shared_ptr<Parameter> para, int level)
{
    this->para = para;
    this->level = level;
    this->kernelUsesFluidNodeIndices = true;

    myPreProcessorTypes.push_back(InitAdvectionDiffusionCompressible);
}

template class F16CompressibleAdvectionDiffusion<TurbulenceModel::None>;
template class F16CompressibleAdvectionDiffusion<TurbulenceModel::Default>;
template class F16CompressibleAdvectionDiffusion<TurbulenceModel::Moeng>;
template class F16CompressibleAdvectionDiffusion<TurbulenceModel::AMDStratified>;

}

//! \}
