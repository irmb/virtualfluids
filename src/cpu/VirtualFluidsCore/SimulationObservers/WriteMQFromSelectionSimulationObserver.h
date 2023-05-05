#ifndef WriteMQFromSelectionSimulationObserver_H
#define WriteMQFromSelectionSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"

#include "LBMSystem.h"
#include "UbTuple.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class WbWriter;
class Block3D;
class GbObject3D;

class WriteMQFromSelectionSimulationObserver : public SimulationObserver
{
public:
    WriteMQFromSelectionSimulationObserver();
    WriteMQFromSelectionSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<GbObject3D> gbObject,
                                    const std::string &path, WbWriter *const writer, SPtr<LBMUnitConverter> conv,
                                    std::shared_ptr<vf::mpi::Communicator> comm);
    ~WriteMQFromSelectionSimulationObserver() override = default;

    void update(real step) override;

protected:
    void collectData(real step);
    void addDataMQ(SPtr<Block3D> block);
    void clearData();

private:
    void init();
    std::vector<UbTupleFloat3> nodes;
    std::vector<std::string> datanames;
    std::vector<std::vector<real>> data;
    std::string path;
    WbWriter *writer;
    SPtr<LBMUnitConverter> conv;
    //   bool bcInformation;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    std::shared_ptr<vf::mpi::Communicator> comm;
    SPtr<GbObject3D> gbObject;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};

#endif
