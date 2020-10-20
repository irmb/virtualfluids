#ifndef WriteGbObjectsCoProcessor_h__
#define WriteGbObjectsCoProcessor_h__

#include "CoProcessor.h"
#include "UbTuple.h"

#include <vector>

class GbObject3D;
class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;

//! \brief     Writes geometry objects as VTK unstructured grid.
//! \details   Writes geometry objects as VTK unstructured grid. Use addGbObject() for add a GbObjects.
//! \author    Konstantin Kutscher
//! \date      December 2018

class WriteGbObjectsCoProcessor : public CoProcessor
{
public:
    WriteGbObjectsCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, WbWriter *const writer,
                              SPtr<Communicator> comm);
    ~WriteGbObjectsCoProcessor() override;
    //! calls collectData.
    void process(double step) override;
    //! adds geometry object
    void addGbObject(SPtr<GbObject3D> object);

protected:
    void collectData(double step);

private:
    std::vector<SPtr<GbObject3D>> objects;
    std::string path;
    WbWriter *writer;
    SPtr<Communicator> comm;
};

#endif // WriteGbObjectsCoProcessor_h__