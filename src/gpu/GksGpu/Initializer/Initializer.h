#ifndef  Initializer_H
#define  Initializer_H

#include <string>
#include <memory>
#include <functional>

#include "VirtualFluidsDefinitions.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "Core/VectorTypes.h"

#include "DataBase/DataBase.h"
#include "FlowStateData/FlowStateData.cuh"

namespace GksGpu {

class VIRTUALFLUIDS_GPU_EXPORT Initializer
{
public:

    static void interpret( SPtr<DataBase> dataBase, std::function<ConservedVariables(Vec3)> initialCondition );

    static void initializeDataUpdate( SPtr<DataBase> dataBase );
};

} // namespace GksGpu

#endif
