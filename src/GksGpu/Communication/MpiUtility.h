#ifndef MpiUtility_H
#define MpiUtility_H

#include <mpi.h>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

class  GksMeshAdapter;

namespace GksGpu {

class  DataBaseAllocator;
struct DataBase;

struct VF_PUBLIC MpiUtility
{
    static int getMpiRankBeforeInit();

    static int getMpiWorldSizeBeforeInit();
};

} // namespace GksGpu

#endif
