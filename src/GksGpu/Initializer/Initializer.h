#ifndef  Initializer_H
#define  Initializer_H

#include <string>
#include <memory>
#include <functional>

#include "VirtualFluidsDefinitions.h"

#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "DataBase/DataBase.h"
#include "FlowStateData/FlowStateData.cuh"

class VF_PUBLIC Initializer
{
public:

    static void interpret( std::shared_ptr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition );

    static Vec3 getCellCenter( std::shared_ptr<DataBase> dataBase, uint cellIdx );

    static void initializeDataUpdate( std::shared_ptr<DataBase> dataBase );
};

#endif
