#include "DistributionDebugWriter.h"
#include "Cuda/CudaMemoryManager.h"
#include "Utilities/testUtilitiesGPU.h"
#include <gmock/gmock.h>

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
