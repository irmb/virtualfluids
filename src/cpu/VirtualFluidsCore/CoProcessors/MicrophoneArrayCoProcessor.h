#ifndef MicrophoneArrayCoProcessor_h__
#define MicrophoneArrayCoProcessor_h__

#include "CoProcessor.h"
#include "LBMSystem.h"
#include "UbTuple.h"
#include <array>
#include <string>
#include <vector>

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class Vector3D;
class DistributionArray3D;

//! \brief     Class implements microphone array.
//! \details   It samples pressure (LBM rho value) and stores to .csv file.
//! \author    Konstantin Kutscher
//! \date      February 2019

class MicrophoneArrayCoProcessor : public CoProcessor
{
public:
    MicrophoneArrayCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                               std::shared_ptr<vf::mpi::Communicator> comm);
    ~MicrophoneArrayCoProcessor() override;

    //! calls collectData.
    void process(real step) override;

    //! add microphone
    bool addMicrophone(Vector3D coords);

protected:
    void collectData(real step);
    void writeFile(real step);

private:
    std::string path;
    std::shared_ptr<vf::mpi::Communicator> comm;

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

#endif // MicrophoneArrayCoProcessor_h__
