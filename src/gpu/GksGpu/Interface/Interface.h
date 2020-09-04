#ifndef  FineToCoarse_H
#define  FineToCoarse_H

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

namespace GksGpu {

class GKSGPU_EXPORT Interface
{
public:
    static void runFineToCoarse( SPtr<DataBase> dataBase, uint level );

    static void runCoarseToFine( SPtr<DataBase> dataBase, uint level );
};

} // namespace GksGpu

#endif
