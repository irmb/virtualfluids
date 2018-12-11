#ifndef  FineToCoarse_H
#define  FineToCoarse_H

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"

class VF_PUBLIC Interface
{
public:
    static void runFineToCoarse( SPtr<DataBase> dataBase, uint level );

    static void runCoarseToFine( SPtr<DataBase> dataBase, uint level );
};

#endif
