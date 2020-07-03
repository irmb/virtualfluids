#ifndef AverageValuesCoProcessor_H
#define AverageValuesCoProcessor_H

#include <PointerDefinitions.h>
#include <vector>
#include <string>

#include "CoProcessor.h"
#include "LBMSystem.h"
#include "UbTuple.h"

class UbScheduler;
class WbWriter;
class Grid3D;
class Block3D;

//! \brief  Computes the time averaged mean velocity and RMS values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s), does averaging according to scheduler (Avs) and resets according to scheduler (rs).  <br>
//!  Computes  the time averaged mean velocity  \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$  and RMS of fluctuations. You need to calculate a square root before plotting RMS. <br>
//           
//! \author  Sonja Uphoff, Kostyantyn Kucher 
// \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$
class AverageValuesCoProcessor : public CoProcessor
{
public:
   AverageValuesCoProcessor();
   AverageValuesCoProcessor(SPtr<Grid3D> grid, const std::string& path, WbWriter* const writer,
       SPtr<UbScheduler> s, SPtr<UbScheduler> Avs, SPtr<UbScheduler> rsMeans, SPtr<UbScheduler> rsRMS, bool restart);
	//! Make update
	void process(double step); 
	//! Resets averaged velocity and RMS-values according to ResetSceduler
	void reset(double step); 
protected:
	//! Prepare data and write in .vtk file
	void collectData(double step);
	//! Reset data
	void resetDataRMS(double step);
	void resetDataMeans(double step);
	//! prepare data
	void addData(const SPtr<Block3D> block);
	void clearData();
	//! Computes average and RMS values of macroscopic quantities 
	void calculateAverageValues(double timeStep);
	////! write .txt file spatial intergrated averaged value, fluctuation, porous features
	//void collectPlotDataZ(double step);
	////! create txt file and write head line 
	//void initPlotDataZ(double step);

private:
	std::vector<UbTupleFloat3> nodes;
	std::vector<UbTupleUInt8> cells;
	std::vector<std::string> datanames;
	std::vector<std::vector<double> > data; 
	std::vector<std::vector<SPtr<Block3D>> > blockVector;
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
	enum Values{AvVx = 0, AvVy = 1, AvVz = 2, AvVxx = 3, AvVyy = 4, AvVzz = 5, AvVxy = 6, AvVxz = 7, AvVyz = 8, AvP = 9, AvPrms = 10}; 

   typedef void (*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/,LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;
};
#endif
