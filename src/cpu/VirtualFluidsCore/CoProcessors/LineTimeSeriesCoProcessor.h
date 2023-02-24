#ifndef LineTimeSeriesCoProcessor_h__
#define LineTimeSeriesCoProcessor_h__

#include <PointerDefinitions.h>
#include <string>

#include <mpi.h>

#include "CoProcessor.h"
#include "LBMSystem.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class GbLine3D;

//! \brief  Writes to .csv file time series for a line in x1 direction.
//! \details It can be used to compute for given time range  the time averaged two-point correlations for a line. <br>
//!  \f$ R_{ij}(x_{a},x{b},t) = <u_{i}(x_{a},t)u_{j}(x_{a}+r,t)> \f$   <br>
//
//! \author  Konstantin Kutscher

class LineTimeSeriesCoProcessor : public CoProcessor
{
public:
    enum Direction { X1, X2, X3 };

public:
    LineTimeSeriesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, SPtr<GbLine3D> line,
                              int level, std::shared_ptr<vf::mpi::Communicator> comm);
    ~LineTimeSeriesCoProcessor() override = default;

    void process(real step) override;
    void writeLine(const std::string &path);

protected:
    void collectData();

private:
    std::string path;
    std::string fname;
    bool root;
    SPtr<GbLine3D> line;
    // function pointer
    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
    int blocknx;
    int blockix1;
    int blockix2;
    int blockix3;
    int level;
    int ix1;
    int ix2;
    int ix3;
    int length;
    MPI_Comm mpi_comm;
    int numOfProc;
    int gridRank;
    Direction dir;
};
#endif // LineTimeSeriesCoProcessor_h__
