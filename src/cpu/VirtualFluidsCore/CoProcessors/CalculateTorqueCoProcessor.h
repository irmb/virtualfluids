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
namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;
class ILBMKernel;

class CalculateTorqueCoProcessor: public CoProcessor 
{
public:
   //! Constructor
   CalculateTorqueCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm);
	virtual ~CalculateTorqueCoProcessor();             
	void process(real step); 
   void addInteractor(SPtr<D3Q27Interactor> interactor);
protected:
	void collectData(real step);
   void calculateForces();
    UbTupleDouble3 getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc);
    UbTupleDouble3 getForcesFromMoments(int x1, int x2, int x3, SPtr<ILBMKernel> kernel, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc, real nx, real ny, real nz);
    UbTupleDouble3 getForcesFromStressTensor(int x1, int x2, int x3, SPtr<ILBMKernel> kernel, SPtr<DistributionArray3D> distributions, SPtr<BoundaryConditions> bc, real nx, real ny, real nz);

private:
   std::string path;
   std::shared_ptr<vf::mpi::Communicator> comm;
   std::vector<SPtr<D3Q27Interactor> > interactors;
   real torqueX1global;
   real torqueX2global;
   real torqueX3global;

   real Fx, Fy, Fz;
};


#endif /* D3Q27ForcesCoProcessor_H */
