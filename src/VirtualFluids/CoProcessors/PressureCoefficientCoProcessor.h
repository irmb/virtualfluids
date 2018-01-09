#ifndef PressureCoefficientCoProcessor_h__
#define PressureCoefficientCoProcessor_h__

#include <memory>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "LBMSystem.h"


class GbCuboid3D;
class D3Q27Interactor;
class Communicator;
class Grid3D;
class UbScheduler;

class PressureCoefficientCoProcessor: public CoProcessor
{
public:
   PressureCoefficientCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s,
       std::shared_ptr<GbCuboid3D> plane, const std::string& path, std::shared_ptr<Communicator> comm);
   ~PressureCoefficientCoProcessor();

   void process(double step) override;

   void addInteractor(std::shared_ptr<D3Q27Interactor> interactor);
   void readValues(int step);

protected:
   void collectData(double step);
   void calculateRho();
   void writeValues(int step);

private:
    std::shared_ptr<GbCuboid3D> plane;
   std::string path;
   std::shared_ptr<Communicator> comm;
   std::vector<std::shared_ptr<D3Q27Interactor> > interactors;
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


