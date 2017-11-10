#ifndef PressureCoefficientCoProcessor_h__
#define PressureCoefficientCoProcessor_h__

#include "CoProcessor.h"
#include "Communicator.h"
#include "GbCuboid3D.h"
#include "D3Q27Interactor.h"

class PressureCoefficientCoProcessor: public CoProcessor
{
public:
   PressureCoefficientCoProcessor(Grid3DPtr grid, UbSchedulerPtr s,
      GbCuboid3DPtr plane, const std::string& path, CommunicatorPtr comm);
   ~PressureCoefficientCoProcessor();
   void process(double step);
   void addInteractor(D3Q27InteractorPtr interactor);
   void readValues(int step);
protected:
   void collectData(double step);
   void calculateRho();
   void writeValues(int step);
private:
   GbCuboid3DPtr plane;
   std::string path;
   CommunicatorPtr comm;
   std::vector<D3Q27InteractorPtr> interactors;
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


