#ifndef Restart_h
#define Restart_h

#include <string>
#include <memory>


#include "GksGpu_export.h"

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"

namespace GksGpu {

struct DataBase;

class GKSGPU_EXPORT Restart
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
