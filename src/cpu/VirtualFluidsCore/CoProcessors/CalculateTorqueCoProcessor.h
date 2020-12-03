/*
 *  D3Q27ForcesCoProcessor.h
 *
 *  Created on: 29.09.2012
 *  Author: K. Kucher
 */

#ifndef CalculateTorqueCoProcessor_H
#define CalculateTorqueCoProcessor_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "CoProcessor.h"
#include "UbTuple.h"
#include "D3Q27System.h"

class ForceCalculator;
class Communicator;
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;

class CalculateTorqueCoProcessor: public CoProcessor 
{
public:
   //! Constructor
   CalculateTorqueCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, SPtr<Communicator> comm);
	virtual ~CalculateTorqueCoProcessor();             
	void process(double step); 
   void addInteractor(SPtr<D3Q27Interactor> interactor);
protected:
	void collectData(double step);
   void calculateForces();
   UbTupleDouble3 getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc);
   void write(std::ofstream *fileObject, double value, char *separator);
private:
   std::string path;
   SPtr<Communicator> comm;
   std::vector<SPtr<D3Q27Interactor> > interactors;
   double forceX1global;
   double forceX2global;
   double forceX3global;

   typedef void(*CalcMacrosFct)(const LBMReal* const& /*f[27]*/, LBMReal& /*rho*/, LBMReal& /*vx1*/, LBMReal& /*vx2*/, LBMReal& /*vx3*/);
   CalcMacrosFct    calcMacrosFct;
   LBMReal compressibleFactor;
};


#endif /* D3Q27ForcesCoProcessor_H */
