#ifndef LineTimeSeriesCoProcessor_h__
#define LineTimeSeriesCoProcessor_h__

#include "CoProcessor.h"
#include "GbLine3D.h"
#include "MPICommunicator.h"

//! \brief  Writes to .csv file time series for a line in x1 direction.
//! \details It can be used to compute for given time range  the time averaged two-point correlations for a line. <br>
//!  \f$ R_{ij}(x_{a},x{b},t) = <u_{i}(x_{a},t)u_{j}(x_{a}+r,t)> \f$   <br>
//           
//! \author  Konstantin Kutscher 

class LineTimeSeriesCoProcessor : public CoProcessor
{
public:
enum Direction {X1, X2, X3};
public:
   LineTimeSeriesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string& path, GbLine3DPtr line, int level, CommunicatorPtr comm);
   ~LineTimeSeriesCoProcessor(){}
   void process(double step);
   void writeLine(const std::string& path);
protected:
   void collectData();
private:
   std::string path;
   std::string fname;
   bool root;
   GbLine3DPtr line;
   //function pointer
   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
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
