#ifndef CellProperties_H
#define CellProperties_H

#ifdef __CUDACC__
#include <cuda_runtime.h>
#else
#define __host__
#define __device__
#endif

//////////////////////////////////////////////////////////////////////////

#define CELL_PROPERTIES_DEFAULT    (0u)
#define CELL_PROPERTIES_GHOST      (1u)
#define CELL_PROPERTIES_WALL       (2u)
#define CELL_PROPERTIES_FINE_GHOST (4u)
#define CELL_PROPERTIES_3          (8u)
#define CELL_PROPERTIES_4          (16u)
#define CELL_PROPERTIES_5          (32u)
#define CELL_PROPERTIES_6          (64u)
#define CELL_PROPERTIES_7          (128u)

//////////////////////////////////////////////////////////////////////////

typedef unsigned char CellProperties;

//////////////////////////////////////////////////////////////////////////

__host__ __device__ inline void setCellProperties( CellProperties& left, const CellProperties& right )
{
    left |= right;
}

__host__ __device__ inline void unsetCellProperties( CellProperties& left, const CellProperties& right )
{
    left &= ~right;
}

__host__ __device__ inline bool isCellProperties( const CellProperties& left, const CellProperties& right )
{
    return (left & right) == right;
}

//////////////////////////////////////////////////////////////////////////

#endif

