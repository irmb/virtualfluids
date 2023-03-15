#ifndef TurbulenceIntensitySimulationObserver_H
#define TurbulenceIntensitySimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "UbTuple.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

class TurbulenceIntensitySimulationObserver : public SimulationObserver
{
public:
    TurbulenceIntensitySimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer,
                                   SPtr<UbScheduler> s, std::shared_ptr<vf::mpi::Communicator> comm);
    void process(real step) override;

protected:
    void collectData(real step);
    void addData(const SPtr<Block3D> block);
    void clearData();
    void calculateAverageValues(real timeStep);

private:
    void init();
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames;
    std::vector<std::vector<real>> data;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::mpi::Communicator> comm;
    enum Values { AvVx = 0, AvVy = 1, AvVz = 2, AvVxxyyzz = 3 };
};
#endif
