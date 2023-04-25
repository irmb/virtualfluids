/*
 *  ForceCalculator.h
 *
 *  Created on: 25.10.2017
 *  Author: S. Peters
 */
#include "lbm/constants/D3Q27.h"
 
#ifndef ForceCalculator_H
#define ForceCalculator_H

#include <memory>
#include <vector>

#include "Vector3D.h"

class D3Q27Interactor;
namespace vf::mpi {class Communicator;}
class DistributionArray3D;
class BoundaryConditions;

class ForceCalculator
{
public:
    ForceCalculator(std::shared_ptr<vf::mpi::Communicator> comm);
    virtual ~ForceCalculator();

    void calculateForces(std::vector<std::shared_ptr<D3Q27Interactor>> interactors);
    Vector3D getForces(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions,
                       std::shared_ptr<BoundaryConditions> bc,
                       const Vector3D &boundaryVelocity = Vector3D(0.0, 0.0, 0.0)) const;

    Vector3D getGlobalForces() const;

private:
    void gatherGlobalForces();

    std::shared_ptr<vf::mpi::Communicator> comm;

    real forceX1global;
    real forceX2global;
    real forceX3global;
};

#endif
