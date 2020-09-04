#ifndef  NestedTimeStep_H
#define  NestedTimeStep_H

#include "VirtualFluidsDefinitions.h"
#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"
namespace GksGpu{ 

class GKSGPU_EXPORT TimeStepping
{
public:

    static void nestedTimeStep( SPtr<DataBase> dataBase, 
                                Parameters parameters, 
                                uint level );

};

} // namespace GksGpu

#endif
