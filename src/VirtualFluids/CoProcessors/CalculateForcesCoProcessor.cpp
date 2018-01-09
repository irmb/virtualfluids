#include "CalculateForcesCoProcessor.h"
#include "BCProcessor.h"

#include "Communicator.h"
#include "D3Q27Interactor.h"
#include "ForceCalculator.h"
#include "UbScheduler.h"
#include "Grid3D.h"

CalculateForcesCoProcessor::CalculateForcesCoProcessor( Grid3DPtr grid, UbSchedulerPtr s, 
                                                    const std::string &path,
                                                    CommunicatorPtr comm ,
                                                    double v, double a, 
                                                    std::shared_ptr<ForceCalculator> forceCalculator) :
                                                    CoProcessor(grid, s),
                                                    path(path), comm(comm),
                                                    v(v), a(a),
                                                    forceCalculator(forceCalculator)
{
   if (comm->getProcessID() == comm->getRoot())
   {
      std::ofstream ostr;
      std::string fname = path;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         std::string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }
      ostr.width(12);
      ostr << "step" << "\t";
      ostr.width(12);
      ostr << "Cx" << "\t";
      ostr.width(12);
      ostr << "Cy"  << "\t"; 
      ostr.width(12);   
      ostr << "Cz" << std::endl;
      ostr.close();
   }
}
//////////////////////////////////////////////////////////////////////////
CalculateForcesCoProcessor::~CalculateForcesCoProcessor()
{

}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::process( double step )
{
   if(scheduler->isDue(step) )
      collectData(step);

   UBLOG(logDEBUG3, "D3Q27ForcesCoProcessor::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void CalculateForcesCoProcessor::collectData( double step )
{
   forceCalculator->calculateForces(this->interactors);

   if (comm->getProcessID() == comm->getRoot())
   {
      int istep = static_cast<int>(step);
      std::ofstream ostr;
      std::string fname = path;
      ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
      if(!ostr)
      { 
         ostr.clear();
         std::string path = UbSystem::getPathFromString(fname);
         if(path.size()>0){ UbSystem::makeDirectory(path); ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);}
         if(!ostr) throw UbException(UB_EXARGS,"couldn't open file "+fname);
      }


      Vector3D globalForces = forceCalculator->getGlobalForces();
      double C1 = getCoefficient(globalForces[0]);
      double C2 = getCoefficient(globalForces[1]);
      double C3 = getCoefficient(globalForces[2]);

      
      ostr.width(12); 
      ostr.setf(std::ios::fixed); 
      ostr << istep << "\t";
      write(&ostr, C1, (char*)"\t");
      write(&ostr, C2, (char*)"\t");
      write(&ostr, C3, (char*)"\t");
      ostr << std::endl;
      ostr.close();
   }
}

double CalculateForcesCoProcessor::getCoefficient(double force)
{
    return 2.0* force / (v*v*a);
}

void CalculateForcesCoProcessor::addInteractor( D3Q27InteractorPtr interactor )
{
   interactors.push_back(interactor);
}

void CalculateForcesCoProcessor::write(std::ofstream *fileObject, double value, char *separator) 
{ 
   (*fileObject).width(12); 
   //(*fileObject).precision(2); 
   (*fileObject).setf(std::ios::fixed); 
   (*fileObject) << value; 
   (*fileObject) << separator; 
} 


