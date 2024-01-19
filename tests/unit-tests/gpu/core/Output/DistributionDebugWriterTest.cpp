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
//! \addtogroup gpu_Output_tests Output
//! \ingroup gpu_core_tests core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include <gmock/gmock.h>

#include <gpu/core/Cuda/CudaMemoryManager.h>
#include <gpu/core/Output/DistributionDebugWriter.h>

#include "../Utilities/testUtilitiesGPU.h"

TEST(DistributionDebugWriterTest, DistributionsAreNotAllocated_CopyDistributions_ShouldThrow)
{
    const auto para = testingVF::createParameterForLevel(0);
    const CudaMemoryManager cudaMemoryManager(para);

    EXPECT_THROW(DistributionDebugWriter::copyDistributionsToHost(*para, cudaMemoryManager), std::runtime_error);
}

TEST(DistributionDebugWriterTest, DistributionsAreNotAllocated_WriteDistributions_ShouldThrow)
{
    const auto para = testingVF::createParameterForLevel(0);
    const CudaMemoryManager cudaMemoryManager(para);

    EXPECT_THROW(DistributionDebugWriter::writeDistributions(*para, 0), std::runtime_error);
}

//! \}
