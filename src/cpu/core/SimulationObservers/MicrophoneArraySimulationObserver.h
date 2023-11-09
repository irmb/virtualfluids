#ifndef MicrophoneArraySimulationObserver_h__
#define MicrophoneArraySimulationObserver_h__

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"
#include <array>
#include <string>
#include <vector>

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class Vector3D;
class DistributionArray3D;

//! \brief     Class implements microphone array.
//! \details   It samples pressure (LBM rho value) and stores to .csv file.
//! \author    Konstantin Kutscher
//! \date      February 2019

class MicrophoneArraySimulationObserver : public SimulationObserver
{
public:
    MicrophoneArraySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                               std::shared_ptr<vf::parallel::Communicator> comm);
    ~MicrophoneArraySimulationObserver() override;

    //! calls collectData.
    void update(real step) override;

    //! add microphone
    bool addMicrophone(Vector3D coords);

protected:
    void collectData(real step);
    void writeFile(real step);

private:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;

    struct Mic {
        unsigned int id;
        SPtr<DistributionArray3D> distridution;
        UbTupleInt3 nodeIndexes;
    };
    std::vector<SPtr<Mic>> microphones;

    std::vector<SPtr<std::stringstream>> strVector;

    int count;
    int micID;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};

#endif // MicrophoneArraySimulationObserver_h__
