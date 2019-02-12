#ifndef  NestedTimeStep_H
#define  NestedTimeStep_H

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

class VF_PUBLIC TimeStepping
{
public:

    static void nestedTimeStep( SPtr<DataBase> dataBase, 
                                Parameters parameters, 
                                uint level );

};

#endif
