#ifndef PressureCoefficientCoProcessor_h__
#define PressureCoefficientCoProcessor_h__

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "LBMSystem.h"
#include "UbTuple.h"


class GbCuboid3D;
class D3Q27Interactor;
class Communicator;
class Grid3D;
class UbScheduler;

class PressureCoefficientCoProcessor: public CoProcessor
{
public:
   PressureCoefficientCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
       SPtr<GbCuboid3D> plane, const std::string& path, SPtr<Communicator> comm);
   ~PressureCoefficientCoProcessor();

   void process(double step) override;

   void addInteractor(SPtr<D3Q27Interactor> interactor);
   void readValues(int step);

protected:
   void collectData(double step);
   void calculateRho();
   void writeValues(int step);

private:
    SPtr<GbCuboid3D> plane;
   std::string path;
   SPtr<Communicator> comm;
   std::vector<SPtr<D3Q27Interactor> > interactors;
   int numberOfSteps;
   double maxStep;

   std::vector<UbTupleFloat3> nodes;
   std::vector<std::string> datanames;
   std::vector<std::vector<double> > data;

   std::vector<double> outValues;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*feq[27]*/, LBMReal& /*(d)rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct calcMacros;
};

#endif // PressureDistributionCoProcessor_h__


