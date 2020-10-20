#ifndef MicrophoneArrayCoProcessor_h__
#define MicrophoneArrayCoProcessor_h__

#include "CoProcessor.h"
#include "LBMSystem.h"
#include "UbTuple.h"
#include <array>
#include <string>
#include <vector>

class Communicator;
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
                               SPtr<Communicator> comm);
    ~MicrophoneArrayCoProcessor() override;

    //! calls collectData.
    void process(double step) override;

    //! add microphone
    bool addMicrophone(Vector3D coords);

protected:
    void collectData(double step);
    void writeFile(double step);

private:
    std::string path;
    SPtr<Communicator> comm;

    struct Mic {
        unsigned int id;
        SPtr<DistributionArray3D> distridution;
        UbTupleInt3 nodeIndexes;
    };
    std::vector<SPtr<Mic>> microphones;

    std::vector<SPtr<std::stringstream>> strVector;

    int count;
    int micID;

    using CalcMacrosFct = void (*)(const LBMReal *const &, LBMReal &, LBMReal &, LBMReal &, LBMReal &);
    CalcMacrosFct calcMacros;
};

#endif // MicrophoneArrayCoProcessor_h__
