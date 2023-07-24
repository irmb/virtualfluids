#ifndef WriteGbObjectsSimulationObserver_h__
#define WriteGbObjectsSimulationObserver_h__

#include "SimulationObserver.h"
#include "UbTuple.h"

#include <vector>

class GbObject3D;
namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;

//! \brief     Writes geometry objects as VTK unstructured grid.
//! \details   Writes geometry objects as VTK unstructured grid. Use addGbObject() for add a GbObjects.
//! \author    Konstantin Kutscher
//! \date      December 2018

class WriteGbObjectsSimulationObserver : public SimulationObserver
{
public:
    WriteGbObjectsSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, WbWriter *const writer,
                              std::shared_ptr<vf::parallel::Communicator> comm);
    ~WriteGbObjectsSimulationObserver() override;
    //! calls collectData.
    void update(real step) override;
    //! adds geometry object
    void addGbObject(SPtr<GbObject3D> object);

protected:
    void collectData(real step);

private:
    std::vector<SPtr<GbObject3D>> objects;
    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::parallel::Communicator> comm;
};

#endif // WriteGbObjectsSimulationObserver_h__