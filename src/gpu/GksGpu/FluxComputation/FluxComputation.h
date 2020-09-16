#ifndef  FluxComputation_H
#define  FluxComputation_H


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

#include "DataBase/DataBase.h"
#include "Parameters/Parameters.h"

namespace GksGpu {

class GKSGPU_EXPORT FluxComputation
{
public:

    static void run( SPtr<DataBase> dataBase, 
                     Parameters parameters, 
                     uint level,
                     bool evaluateCommFaces = false);
};

} // namespace GksGpu

#endif
