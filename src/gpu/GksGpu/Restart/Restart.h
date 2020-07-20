#ifndef Restart_h
#define Restart_h

#include <string>
#include <memory>

#include <VirtualFluidsDefinitions.h>

#include "Core/PointerDefinitions.h"
#include "Core/DataTypes.h"

namespace GksGpu {

struct DataBase;

class VIRTUALFLUIDS_GPU_EXPORT Restart
{

public:
    static void writeRestart( SPtr<DataBase> dataBase, std::string filename, uint  iter );

    static bool readRestart ( SPtr<DataBase> dataBase, std::string filename, uint& iter );

private:
    Restart(){}
    ~Restart(){}
};

} // namespace GksGpu


#endif
