#ifndef D3Q27MEANVALUESPOSTPROCESSOR_H
#define D3Q27MEANVALUESPOSTPROCESSOR_H

#include <Postprocessor.h>
#include <D3Q27IntegrateValuesHelper.h>
#include <LBMUnitConverter.h>

class D3Q27MeanValuesPostprocessor : public Postprocessor
{
public:
   D3Q27MeanValuesPostprocessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string &path,
      CommunicatorPtr comm, D3Q27IntegrateValuesHelperPtr ih, double factorPressure, double factorVelocity);
                                
   virtual ~D3Q27MeanValuesPostprocessor();
   void update(double step);
   
private:
   D3Q27IntegrateValuesHelperPtr ih;
   CommunicatorPtr comm;
   std::string filename;
   double factorPressure, factorVelocity;
   void collectPostprocessData(double step);
};

#endif
   
