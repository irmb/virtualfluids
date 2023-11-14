#ifndef DISTRIBUTION_DEBUG_WRITER
#define DISTRIBUTION_DEBUG_WRITER

#include <basics/DataTypes.h>

class Parameter;
class CudaMemoryManager;

//! \brief Functions to write the distributions ("f's") to a VTK-file
//! \details to make this work, the distributions need to be copied to the host. This can be slow!
class DistributionDebugWriter
{
public:
    //! \brief allocate memory for the distributions on the host
    //! \details Call only once. If no memory is allocated, error
    static void allocateDistributionsOnHost(const CudaMemoryManager& cudaMemoryManager);
    //! \param level allocate memory for the specified level only
    static void allocateDistributionsOnHost(const CudaMemoryManager& cudaMemoryManager, uint level);

    //! \brief Copy distributions from device to host. Call this function before writing data
    //! \details copies data for all levels
    static void copyDistributionsToHost(const Parameter& para, const CudaMemoryManager& cudaMemoryManager);
    //! \param level allocate memory for the specified level only
    static void copyDistributionsToHost(const Parameter& para, const CudaMemoryManager& cudaMemoryManager, uint level);

    //! \brief write the distributions for all levels
    //! \details for this to work the distributions have to be allocated on the host and have to be copied to the host for
    //! the current timestep
    static void writeDistributions(const Parameter& para, uint timestep);
    //! \brief write the distributions for one level
    //! \details for this to work the distributions have to be allocated on the host and have to be copied to the host for
    //! the current timestep
    static void writeDistributionsForLevel(const Parameter& para, uint level, uint timestep);
};

#endif
