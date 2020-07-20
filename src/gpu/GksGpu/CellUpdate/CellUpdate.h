#ifndef  CellUpdate_H
#define  CellUpdate_H

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

namespace GksGpu {

class VIRTUALFLUIDS_GPU_EXPORT CellUpdate
{
public:

    static void run( SPtr<DataBase> dataBase, 
                     Parameters parameters, 
                     uint level );
};

} // namespace GksGpu

#endif
