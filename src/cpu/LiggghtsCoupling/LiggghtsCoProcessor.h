//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LiggghtsCoProcessor.h
//! \ingroup LiggghtsCoupling
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LIGGGHTS_CO_PROCESSOR_H
#define LIGGGHTS_CO_PROCESSOR_H

#include "CoProcessor.h"

#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

class LiggghtsCoProcessor : public CoProcessor
{
public:
    LiggghtsCoProcessor();
    virtual ~LiggghtsCoProcessor();

};

#endif

///*
// *  Author: S. Peters
// *  mail: peters@irmb.tu-bs.de
// */
//#ifndef DEM_CO_PROCESSOR_H
//#define DEM_CO_PROCESSOR_H
//
//#include <map>
//#include <memory>
//#include <vector>
//
//#include "Vector3D.h"
//
//#include "CoProcessor.h"
//#include "UbTuple.h"
//
//#include <pe/basic.h>
//
////#define TIMING
//
//#ifdef TIMING
//#include "UbTiming.h"
//#endif
//
//class PhysicsEngineGeometryAdapter;
//class PhysicsEngineSolverAdapter;
//class PePhysicsEngineSolverAdapter;
//class PhysicsEngineMaterialAdapter;
//class PePhysicsEngineGeometryAdapter;
//
//class UbScheduler;
//class Grid3D;
//class ForceCalculator;
//class Communicator;
//class MovableObjectInteractor;
//class Communicator;
//class BoundaryConditionsBlockVisitor;
//
//class DemCoProcessor : public CoProcessor
//{
//public:
//    DemCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s, std::shared_ptr<Communicator> comm,
//                   std::shared_ptr<ForceCalculator> forceCalculator,
//                   std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver, double intermediatePeSteps = 1.0);
//    virtual ~DemCoProcessor();
//
//    void addInteractor(std::shared_ptr<MovableObjectInteractor> interactor,
//                       std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial,
//                       Vector3D initalVelocity = Vector3D(0.0, 0.0, 0.0));
//    void process(double step) override;
//    std::shared_ptr<PhysicsEngineSolverAdapter> getPhysicsEngineSolver();
//    void distributeIDs();
//    void setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> blockVisitor);
//    bool isDemObjectInAABB(std::array<double, 6> AABB);
//    int addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles);
//    void getObjectsPropertiesVector(std::vector<double> &p);
//    void addPeGeo(walberla::pe::RigidBody *peGeo);
//    void removePeGeo(walberla::pe::RigidBody *peGeo);
//    void addPeShadowGeo(walberla::pe::RigidBody *peGeo);
//    void removePeShadowGeo(walberla::pe::RigidBody *peGeo);
//    bool isSpheresIntersection(double centerX1, double centerX2, double centerX3, double d);
//
//private:
//    std::shared_ptr<PhysicsEngineGeometryAdapter>
//    createPhysicsEngineGeometryAdapter(std::shared_ptr<MovableObjectInteractor> interactor,
//                                       std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const;
//    void applyForcesOnGeometries();
//    void setForcesToObject(SPtr<Grid3D> grid, std::shared_ptr<MovableObjectInteractor> interactor,
//                           std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry);
//    void scaleForcesAndTorques(double scalingFactor);
//    void calculateDemTimeStep(double step);
//    void moveVfGeoObjects();
//    walberla::pe::RigidBody *getPeGeoObject(walberla::id_t id);
//    std::shared_ptr<PePhysicsEngineGeometryAdapter> getPeGeoAdapter(unsigned long long systemId);
//
//private:
//    std::shared_ptr<Communicator> comm;
//    std::vector<std::shared_ptr<MovableObjectInteractor>> interactors;
//    std::shared_ptr<ForceCalculator> forceCalculator;
//    std::shared_ptr<PePhysicsEngineSolverAdapter> physicsEngineSolver;
//    std::vector<std::shared_ptr<PhysicsEngineGeometryAdapter>> physicsEngineGeometrieAdapters;
//    double intermediateDemSteps;
//    SPtr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor;
//    // walberla::pe::BodyStorage* bodyStorage;    //!< Reference to the central body storage.
//    // walberla::pe::BodyStorage* bodyStorageShadowCopies;    //!< Reference to the body storage containing body shadow
//    // copies.
//
//    std::map<unsigned long long, std::shared_ptr<PePhysicsEngineGeometryAdapter>> geoIdMap;
//
//#ifdef TIMING
//    UbTimer timer;
//#endif
//};
//
//#endif
