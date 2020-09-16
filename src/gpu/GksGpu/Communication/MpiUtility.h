#ifndef MpiUtility_H
#define MpiUtility_H

#include <mpi.h>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

class  GksMeshAdapter;

namespace GksGpu {

class  DataBaseAllocator;
struct DataBase;

struct GKSGPU_EXPORT MpiUtility
{
    static int getMpiRankBeforeInit();

    static int getMpiWorldSizeBeforeInit();
};

} // namespace GksGpu

#endif
