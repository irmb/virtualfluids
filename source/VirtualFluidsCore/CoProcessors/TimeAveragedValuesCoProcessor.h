#ifndef TimeAveragedValuesCoProcessor_H
#define TimeAveragedValuesCoProcessor_H

#include "CoProcessor.h"
#include "Grid3D.h"
#include "Block3D.h"
#include "LBMUnitConverter.h"
#include "Communicator.h"
#include "D3Q27IntegrateValuesHelper.h"
#include "WbWriter.h"

#include <boost/shared_ptr.hpp>
class TimeAveragedValuesCoProcessor;
typedef boost::shared_ptr<TimeAveragedValuesCoProcessor> TimeAveragedValuesCoProcessorPtr;

//! \brief  Computes the time averaged mean velocity and RMS values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s), does averaging according to scheduler (Avs) and resets according to scheduler (rs).  <br>
//!  Computes  the time averaged mean velocity  \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$  and RMS of fluctuations. You need to calculate a square root before plotting RMS. <br>
//           
//! \author  Konstantin Kutscher 
// \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$

//struct plotZ
//{
//
//};

class TimeAveragedValuesCoProcessor : public CoProcessor
{
public:
   enum Options
   {
      Velocity = 1,
      Fluctuations = 2,
      Triplecorrelations = 4
   };
public:
   TimeAveragedValuesCoProcessor();
   TimeAveragedValuesCoProcessor(Grid3DPtr grid, const std::string& path, WbWriter* const writer,
      UbSchedulerPtr s, int options);
   //! Make update
   void process(double step);
   //! Resets averaged velocity and RMS-values according to ResetSceduler
   void reset(double step);
protected:
   //! Prepare data and write in .vtk file
   void collectData(double step);
   //! prepare data
   void addData(const Block3DPtr block);
   void clearData();
   //! Computes average values of velocity , fluctuations and triple correlations 
   void calculateAverageValues(double timeStep);
   //! Computes subtotal of velocity , fluctuations and triple correlations
   void calculateSubtotal();

private:
   std::vector<UbTupleFloat3> nodes;
   std::vector<UbTupleInt8> cells;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;
   std::vector<std::vector<Block3DPtr> > blockVector;
   int minInitLevel; //min init level
   int maxInitLevel;
   int gridRank;
   int resetStepRMS;
   int resetStepMeans;
   double averageInterval;
   std::string path;
   WbWriter* writer;
   bool restart, compressible;
   UbSchedulerPtr averageScheduler;  //additional scheduler to averaging after a given interval
   UbSchedulerPtr resetSchedulerRMS;  //additional scheduler to restart averaging after a given interval
   UbSchedulerPtr resetSchedulerMeans;  //additional scheduler to restart averaging after a given interval
   //labels for the different components, e.g. AvVxx for time averaged RMS: 1/n SUM((U-Umean)^2)
   //you need to calculate a square root before plotting RMS
   enum Velocity { Vx, Vy, Vz };
   enum Fluctuations { Vxx, Vyy, Vzz, Vxy, Vxz, Vyz };
   enum Triplecorrelations { Vxxx, Vxxy, Vxxz, Vyyy, Vyyx, Vyyz, Vzzz, Vzzx, Vzzy, Vxyz };
   //enum Pressure { P, Prms };

   int options;
   int counter;
   double breakStep;

   int iMinX1, iMinX2, iMinX3;
   int iMaxX1, iMaxX2, iMaxX3;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;

   

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
