#include "D3Q27MeanValuesPostprocessor.h"

D3Q27MeanValuesPostprocessor::D3Q27MeanValuesPostprocessor(Grid3DPtr grid, UbSchedulerPtr s, const std::string &path,
                                                           CommunicatorPtr comm, D3Q27IntegrateValuesHelperPtr ih, double factorPressure, double factorVelocity) :
                                                            Postprocessor(grid, s),
                                                            filename(path),
                                                            comm(comm),
                                                            ih(ih),
                                                            factorPressure(factorPressure),
                                                            factorVelocity(factorVelocity)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      FILE *output;
      output = fopen(filename.c_str(), "w");
      if (output == NULL)
      {
         std::string pathf = UbSystem::getPathFromString(filename);
         if (filename.size()>0) { UbSystem::makeDirectory(pathf); output = fopen(filename.c_str(), "w"); }
         if (output == NULL) throw UbException(UB_EXARGS, "can not open " + filename);
      }
         
      //fprintf(output, "%-18s %-18s %-18s %-18s %-18s %-18s %-18s\n", "time step", "area", "pressure [Pa*m^2]", "velocity [m^3/s]", "velocityX [m^3/s]", "velocityY [m^3/s]", "velocityZ [m^3/s]");
      fprintf(output, "%-18s;%-18s;%-18s;%-18s;%-18s\n", "time step", "area", "velocityX [m^3/s]", "velocityY [m^3/s]", "velocityZ [m^3/s]");
      fclose(output);
   }
}
//////////////////////////////////////////////////////////////////////////
D3Q27MeanValuesPostprocessor::~D3Q27MeanValuesPostprocessor()
{

}
//////////////////////////////////////////////////////////////////////////
void D3Q27MeanValuesPostprocessor::update(double step)
{
   if (scheduler->isDue(step))
      collectPostprocessData(step);
}
//////////////////////////////////////////////////////////////////////////
void D3Q27MeanValuesPostprocessor::collectPostprocessData(double step)
{
   ih->calculateMQ();
   
   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = static_cast<int>(step);
      std::ofstream ostr;
      double nn = ih->getNumberOfFluidsNodes();
      //double p = (ih->getPress())*factorPressure;
      double vx1 = (ih->getVx1())*factorVelocity;
      double vx2 = (ih->getVx2())*factorVelocity;
      double vx3 = (ih->getVx3())*factorVelocity;
      //double vm = (ih->getVm())*factorVelocity;

      FILE *output;
      output = fopen(filename.c_str(), "a");
      if (output == NULL)
         throw UbException(UB_EXARGS, "can not open " + filename);

      //fprintf(output, "%-18d %-18d %-18.3f %-18.3f %-18.3f %-18.3f %-18.3f\n", istep, (int)nn, p, vm, vx1, vx2, vx3);
      fprintf(output, "%-18d %-18d %-18.3f %-18.3f %-18.3f\n", istep, (int)nn, vx1, vx2, vx3);

      fclose(output);

      UBLOG(logINFO, "D3Q27MeanValuesPostprocessor step: " << istep);
   }
}
