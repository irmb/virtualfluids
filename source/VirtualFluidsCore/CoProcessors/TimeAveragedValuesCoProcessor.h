#ifndef TimeAveragedValuesCoProcessor_H
#define TimeAveragedValuesCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "LBMSystem.h"

class Communicator;
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

//! \brief  Computes the time averaged mean velocity and RMS values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s), does averaging according to scheduler (Avs) and resets according to scheduler (rs).  <br>
//!  Computes  the time averaged mean velocity  \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$  and RMS of fluctuations. You need to calculate a square root before plotting RMS. <br>
//           
//! \author  Konstantin Kutscher 
// \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$


class TimeAveragedValuesCoProcessor : public CoProcessor
{
public:
   enum Options
   {
      Density            = 1,
      Velocity           = 2,
      Fluctuations       = 4,
      Triplecorrelations = 8,

      //Velocity           = 1,
      //Fluctuations       = 2,
      //Triplecorrelations = 4,
   };
public:
   TimeAveragedValuesCoProcessor();
   TimeAveragedValuesCoProcessor(SPtr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       SPtr<UbScheduler> s, SPtr<Communicator> comm, int options);
   TimeAveragedValuesCoProcessor(SPtr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       SPtr<UbScheduler> s, SPtr<Communicator> comm, int options, std::vector<int> levels, std::vector<double>& levelCoords, std::vector<double>& bounds, bool timeAveraging = true);
   //! Make update
   void process(double step);
   //! Computes subtotal of velocity , fluctuations and triple correlations
   void calculateSubtotal(double step);
   void addLevelCoordinate(double c);
   void reset();
   void setWithGhostLayer(bool val);
   bool getWithGhostLayer();

protected:
   //! Prepare data and write in .vtk file
   void collectData(double step);
   //! prepare data
   void addData(const SPtr<Block3D> block);
   void clearData();
   //! Computes average values of velocity , fluctuations and triple correlations 
   void calculateAverageValues(double timeStep);

   void init(SPtr<UbScheduler> s);
   void planarAverage(double step);

private:
    SPtr<Communicator> comm;
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::vector<std::vector<SPtr<Block3D> > > blockVector;
   bool root;
   int minInitLevel; //min init level
   int maxInitLevel;
   int gridRank;
   int resetStepRMS;
   int resetStepMeans;
   double averageInterval;
   std::string path;
   WbWriter* writer;
   bool restart, compressible;
   SPtr<UbScheduler> averageScheduler;  //additional scheduler to averaging after a given interval
   SPtr<UbScheduler> resetSchedulerRMS;  //additional scheduler to restart averaging after a given interval
   SPtr<UbScheduler> resetSchedulerMeans;  //additional scheduler to restart averaging after a given interval
   //labels for the different components, e.g. AvVxx for time averaged RMS: 1/n SUM((U-Umean)^2)
   //you need to calculate a square root before plotting RMS
   enum Velocity { Vx, Vy, Vz };
   enum Fluctuations { Vxx, Vyy, Vzz, Vxy, Vxz, Vyz };
   enum Triplecorrelations { Vxxx, Vxxy, Vxxz, Vyyy, Vyyx, Vyyz, Vzzz, Vzzx, Vzzy, Vxyz };
   //enum Pressure { P, Prms };

   int options;
   int lcounter;
   int numberOfFineSteps;
   int numberOfSteps;
   int fineStep;
   int minFineStep;
   int maxFineStep;
   int levelFactor;
   double maxStep;

   int iMinX1, iMinX2, iMinX3;
   //int iMaxX1, iMaxX2, iMaxX3;
   int iMinC;
   int iMaxC;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;

   bool planarAveraging;
   bool timeAveraging;
   std::vector<double> levelCoords;
   std::vector<int> levels;
   std::vector<double> bounds;

   bool withGhostLayer;

   //friend class boost::serialization::access;
   //template<class Archive>
   //void serialize(Archive & ar, const unsigned int version)
   //{
   //   ar & boost::serialization::base_object<CoProcessor>(*this);
   //   ar & path;
   //   ar & blockVector;
   //   ar & minInitLevel;
   //   ar & maxInitLevel;
   //   ar & gridRank;
   //   ar & writer;
   //   ar & resetStepRMS;
   //   ar & resetStepMeans;
   //   ar & AverageInterval;
   //   ar & averageScheduler;
   //   ar & resetSchedulerRMS;
   //   ar & resetSchedulerMeans;
   //   ar & compressible;
   //}
};
#endif
