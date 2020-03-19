#ifndef  Initializer_H
#define  Initializer_H

#include <string>
#include <memory>
#include <functional>

#include "VirtualFluidsDefinitions.h"

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "DataBase/DataBase.h"
#include "FlowStateData/FlowStateData.cuh"

namespace GksGpu {

class VF_PUBLIC Initializer
{
public:

    static void interpret( SPtr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition );

    static void initializeDataUpdate( SPtr<DataBase> dataBase );
};

} // namespace GksGpu

#endif
