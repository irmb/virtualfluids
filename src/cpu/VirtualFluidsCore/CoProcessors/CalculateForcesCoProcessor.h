/*
 *  D3Q27ForcesCoProcessor.h
 *
 *  Created on: 29.09.2012
 *  Author: K. Kucher
 */

#ifndef D3Q27ForcesCoProcessor_H
#define D3Q27ForcesCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "UbTuple.h"

class ForceCalculator;
class Communicator;
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;

class CalculateForcesCoProcessor: public CoProcessor 
{
public:
   //! Constructor
   //! \param v - velocity of fluid in LB units
   //! \param a - area of object in LB units
   CalculateForcesCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
       SPtr<Communicator> comm, double v, double a);
	~CalculateForcesCoProcessor() override;             
	void process(double step) override; 
   void addInteractor(SPtr<D3Q27Interactor> interactor);
protected:
	void collectData(double step);
   void calculateForces();
   UbTupleDouble3 getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc);
   void calculateCoefficients();
   void write(std::ofstream *fileObject, double value, char *separator);
private:
   std::string path;
   SPtr<Communicator> comm;
   std::vector<SPtr<D3Q27Interactor> > interactors;
   double forceX1global;
   double forceX2global;
   double forceX3global;
   double v;     //!< is the speed of the object relative to the fluid
   double a;     //!< is the reference area
   double C1;
   double C2;
   double C3;
};


#endif /* D3Q27ForcesCoProcessor_H */
