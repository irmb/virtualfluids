#ifndef D3Q27MEANVALUESPOSTPROCESSOR_H
#define D3Q27MEANVALUESPOSTPROCESSOR_H

#include <CoProcessor.h>
#include <D3Q27IntegrateValuesHelper.h>
#include <LBMUnitConverter.h>

class MeanValuesCoProcessor : public CoProcessor
{
public:
   MeanValuesCoProcessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string &path,
      CommunicatorPtr comm, D3Q27IntegrateValuesHelperPtr ih, double factorPressure, double factorVelocity);
                                
   virtual ~MeanValuesCoProcessor();
   void process(double step);
   
private:
   D3Q27IntegrateValuesHelperPtr ih;
   CommunicatorPtr comm;
   std::string filename;
   double factorPressure, factorVelocity;
   void collectPostprocessData(double step);
};

#endif
   
