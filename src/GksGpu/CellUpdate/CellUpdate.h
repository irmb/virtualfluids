#ifndef  CellUpdate_H
#define  CellUpdate_H

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

class VF_PUBLIC CellUpdate
{
public:

    static void updateCells( SPtr<DataBase> dataBase, 
                             Parameters parameters, 
                             uint level );
};

#endif
