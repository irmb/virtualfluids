#include "LiggghtsCoProcessor.h"
//
//#include "BCProcessor.h"
//#include "Communicator.h"
//#include "DataSet3D.h"
//#include "DistributionArray3D.h"
//#include "ForceCalculator.h"
//#include "GbSphere3D.h"
//#include "Grid3D.h"
//#include "ILBMKernel.h"
//#include "MovableObjectInteractor.h"
//#include "SetBcBlocksBlockVisitor.h"
//#include "UbScheduler.h"
//
//#include "PePhysicsEngineGeometryAdapter.h"
//#include "PePhysicsEngineSolverAdapter.h"
//#include "PhysicsEngineGeometryAdapter.h"
//#include "PhysicsEngineMaterialAdapter.h"
//#include "PhysicsEngineSolverAdapter.h"
//
//#include "BCArray3D.h"
//#include "Block3D.h"
//#include "BoundaryConditions.h"
//#include "BoundaryConditionsBlockVisitor.h"
//#include "MPICommunicator.h"
//
//#include "UbLogger.h"
//
//#include <array>
//#include <functional>
//
//DemCoProcessor::DemCoProcessor(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<Communicator> comm,
//                               std::shared_ptr<ForceCalculator> forceCalculator,
//                               std::shared_ptr<PhysicsEngineSolverAdapter> physicsEngineSolver,
//                               double intermediatePeSteps)
//    : CoProcessor(grid, s), comm(comm), forceCalculator(forceCalculator),
//      physicsEngineSolver(std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)),
//      intermediateDemSteps(intermediatePeSteps)
//{
//#ifdef TIMING
//    timer.resetAndStart();
//#endif
//
//    std::shared_ptr<walberla::blockforest::BlockForest> forest =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
//    std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();
//
//    for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt) {
//        walberla::pe::Storage *storage                     = blockIt->getData<walberla::pe::Storage>(*storageId.get());
//        walberla::pe::BodyStorage *bodyStorage             = &(*storage)[0];
//        walberla::pe::BodyStorage *bodyStorageShadowCopies = &(*storage)[1];
//
//        bodyStorage->registerAddCallback("DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::addPeGeo), this));
//        bodyStorage->registerRemoveCallback("DemCoProcessor",
//                                            std::bind1st(std::mem_fun(&DemCoProcessor::removePeGeo), this));
//
//        bodyStorageShadowCopies->registerAddCallback("DemCoProcessor",
//                                                     std::bind1st(std::mem_fun(&DemCoProcessor::addPeShadowGeo), this));
//        bodyStorageShadowCopies->registerRemoveCallback(
//            "DemCoProcessor", std::bind1st(std::mem_fun(&DemCoProcessor::removePeShadowGeo), this));
//    }
//}
//
//DemCoProcessor::~DemCoProcessor()
//{
//    std::shared_ptr<walberla::blockforest::BlockForest> forest =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
//    std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();
//
//    for (auto &currentBlock : *forest) {
//        walberla::pe::Storage *storage           = currentBlock.getData<walberla::pe::Storage>(*storageId.get());
//        walberla::pe::BodyStorage &localStorage  = (*storage)[0];
//        walberla::pe::BodyStorage &shadowStorage = (*storage)[1];
//
//        localStorage.clearAddCallbacks();
//        localStorage.clearRemoveCallbacks();
//
//        shadowStorage.clearAddCallbacks();
//        shadowStorage.clearRemoveCallbacks();
//    }
//}
//
//void DemCoProcessor::addInteractor(std::shared_ptr<MovableObjectInteractor> interactor,
//                                   std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial,
//                                   Vector3D initalVelocity)
//{
//    interactors.push_back(interactor);
//    const int id = static_cast<int>(interactors.size() - 1);
//    interactor->setID(id);
//    const auto peGeometryAdapter = this->createPhysicsEngineGeometryAdapter(interactor, physicsEngineMaterial);
//    if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometryAdapter)->isActive()) {
//        peGeometryAdapter->setLinearVelolocity(initalVelocity);
//        geoIdMap.insert(
//            std::make_pair(std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometryAdapter)->getSystemID(),
//                           std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(peGeometryAdapter)));
//    }
//    SetBcBlocksBlockVisitor setBcVisitor(interactor);
//    grid->accept(setBcVisitor);
//
//    // std::vector< std::shared_ptr<Block3D> > blockVector;
//    // UbTupleInt3 blockNX=grid->getBlockNX();
//    // SPtr<GbObject3D> geoObject(interactor->getGbObject3D());
//    // double ext = 0.0;
//    // std::array<double, 6> AABB ={
//    // geoObject->getX1Minimum(),geoObject->getX2Minimum(),geoObject->getX3Minimum(),geoObject->getX1Maximum(),geoObject->getX2Maximum(),geoObject->getX3Maximum()
//    // }; grid->getBlocksByCuboid(AABB[0]-(double)val<1>(blockNX)*ext, AABB[1]-(double)val<2>(blockNX)*ext,
//    // AABB[2]-(double)val<3>(blockNX)*ext, AABB[3]+(double)val<1>(blockNX)*ext, AABB[4]+(double)val<2>(blockNX)*ext,
//    // AABB[5]+(double)val<3>(blockNX)*ext, blockVector); for (std::shared_ptr<Block3D> block : blockVector)
//    //{
//    //   if (block->getKernel())
//    //   {
//    //      interactor->setBCBlock(block);
//    //      //UBLOG(logINFO, "DemCoProcessor::addInteractor() rank = "<<comm->getProcessID());
//    //   }
//    //}
//
//    interactor->initInteractor();
//
//    physicsEngineGeometrieAdapters.push_back(peGeometryAdapter);
//}
//
//std::shared_ptr<PhysicsEngineGeometryAdapter> DemCoProcessor::createPhysicsEngineGeometryAdapter(
//    std::shared_ptr<MovableObjectInteractor> interactor,
//    std::shared_ptr<PhysicsEngineMaterialAdapter> physicsEngineMaterial) const
//{
//    const int id              = static_cast<int>(interactors.size() - 1);
//    SPtr<GbSphere3D> vfSphere = std::static_pointer_cast<GbSphere3D>(interactor->getGbObject3D());
//    const Vector3D position(vfSphere->getX1Centroid(), vfSphere->getX2Centroid(), vfSphere->getX3Centroid());
//    auto peGeometryAdapter = this->physicsEngineSolver->createPhysicsEngineGeometryAdapter(
//        id, position, vfSphere->getRadius(), physicsEngineMaterial);
//    interactor->setPhysicsEngineGeometry(peGeometryAdapter);
//    return peGeometryAdapter;
//}
//
//void DemCoProcessor::process(double actualTimeStep)
//{
//#ifdef TIMING
//    timer.resetAndStart();
//#endif
//
//    this->applyForcesOnGeometries();
//
//#ifdef TIMING
//    if (comm->isRoot())
//        UBLOG(logINFO, "DemCoProcessor::process start step: " << actualTimeStep);
//    if (comm->isRoot())
//        UBLOG(logINFO, "DemCoProcessor::applyForcesOnGeometries() time = " << timer.stop() << " s");
//#endif
//
//    if (scheduler->isDue(actualTimeStep)) {
//        // UBLOG(logINFO, "DemCoProcessor::update - START - timestep = " << actualTimeStep);
//        const double demTimeStepsPerIteration = scheduler->getMinStep();
//
//        if (demTimeStepsPerIteration != 1)
//            this->scaleForcesAndTorques(1.0 / demTimeStepsPerIteration);
//
//#ifdef TIMING
//        if (comm->isRoot())
//            UBLOG(logINFO, "DemCoProcessor::scaleForcesAndTorques() time = " << timer.stop() << " s");
//        if (comm->isRoot())
//            UBLOG(logINFO, "DemCoProcessor::calculateDemTimeStep():");
//#endif
//
//        if (this->intermediateDemSteps == 1)
//            this->calculateDemTimeStep(demTimeStepsPerIteration);
//
//        //#ifdef TIMING
//        //      if (comm->isRoot()) UBLOG(logINFO, "DemCoProcessor::calculateDemTimeStep() time = "<<timer.stop()<<"
//        //      s");
//        //#endif
//        // if ((int)actualTimeStep % 100 == 0)
//        //{
//        //    if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometries[0])->isActive())
//        //    {
//        //        //UBLOG(logINFO, "v: (x,y,z) " << physicsEngineGeometries[0]->getLinearVelocity() << " actualTimeStep
//        //        = " << UbSystem::toString(actualTimeStep));
//        //    }
//        //}
//
//        // during the intermediate time steps of the collision response, the currently acting forces
//        // (interaction forces, gravitational force, ...) have to remain constant.
//        // Since they are reset after the call to collision response, they have to be stored explicitly before.
//        // Then they are set again after each intermediate step.
//
//        this->moveVfGeoObjects();
//
//#ifdef TIMING
//        if (comm->isRoot())
//            UBLOG(logINFO, "DemCoProcessor::moveVfGeoObject() time = " << timer.stop() << " s");
//#endif
//
//        grid->accept(*boundaryConditionsBlockVisitor.get());
//
//#ifdef TIMING
//        if (comm->isRoot())
//            UBLOG(logINFO, "grid->accept(*boundaryConditionsBlockVisitor.get()) time = " << timer.stop() << " s");
//#endif
//
//        // UBLOG(logINFO, "DemCoProcessor::update - END - timestep = " << actualTimeStep);
//    }
//
//#ifdef TIMING
//    if (comm->isRoot())
//        UBLOG(logINFO, "DemCoProcessor::process stop step: " << actualTimeStep);
//#endif
//}
////////////////////////////////////////////////////////////////////////////
//std::shared_ptr<PhysicsEngineSolverAdapter> DemCoProcessor::getPhysicsEngineSolver() { return physicsEngineSolver; }
//
//void DemCoProcessor::applyForcesOnGeometries()
//{
//    for (int i = 0; i < physicsEngineGeometrieAdapters.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            this->setForcesToObject(grid, interactors[i], physicsEngineGeometrieAdapters[i]);
//
//            // physicsEngineGeometries[i]->setLinearVelolocity(Vector3D(-0.001, 0.0, 0.0));
//            // physicsEngineGeometries[i]->setAngularVelocity(Vector3D(0.01, 0.01, 0.01));
//            // UBLOG(logINFO, "v: (x,y,z) " << physicsEngineGeometries[i]->getLinearVelocity());
//        }
//    }
//}
//
//void DemCoProcessor::setForcesToObject(SPtr<Grid3D> grid, SPtr<MovableObjectInteractor> interactor,
//                                       std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry)
//{
//    for (BcNodeIndicesMap::value_type t : interactor->getBcNodeIndicesMap()) {
//        SPtr<Block3D> block                     = t.first;
//        SPtr<ILBMKernel> kernel                 = block->getKernel();
//        SPtr<BCArray3D> bcArray                 = kernel->getBCProcessor()->getBCArray();
//        SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
//        distributions->swap();
//
//        std::set<std::vector<int>> &transNodeIndicesSet = t.second;
//        for (std::vector<int> node : transNodeIndicesSet) {
//            int x1 = node[0];
//            int x2 = node[1];
//            int x3 = node[2];
//
//            if (kernel->isInsideOfDomain(x1, x2, x3) && bcArray->isFluid(x1, x2, x3)) {
//                // TODO: calculate assumed boundary position
//
//                const Vector3D worldCoordinates = grid->getNodeCoordinates(block, x1, x2, x3);
//                const auto boundaryVelocity     = physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);
//
//                SPtr<BoundaryConditions> bc = bcArray->getBC(x1, x2, x3);
//                const Vector3D force = forceCalculator->getForces(x1, x2, x3, distributions, bc, boundaryVelocity);
//                physicsEngineGeometry->addForceAtPosition(force, worldCoordinates);
//            }
//        }
//        distributions->swap();
//    }
//}
//
//void DemCoProcessor::scaleForcesAndTorques(double scalingFactor)
//{
//    for (int i = 0; i < physicsEngineGeometrieAdapters.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            const Vector3D force  = physicsEngineGeometrieAdapters[i]->getForce() * scalingFactor;
//            const Vector3D torque = physicsEngineGeometrieAdapters[i]->getTorque() * scalingFactor;
//
//            physicsEngineGeometrieAdapters[i]->resetForceAndTorque();
//
//            physicsEngineGeometrieAdapters[i]->setForce(force);
//            physicsEngineGeometrieAdapters[i]->setTorque(torque);
//
//            // UBLOG(logINFO, "F: (x,y,z) " << force);
//            // UBLOG(logINFO, "T: (x,y,z) " << torque);
//        }
//    }
//}
//
//void DemCoProcessor::calculateDemTimeStep(double step)
//{
//    physicsEngineSolver->runTimestep(step);
//
//#ifdef TIMING
//    if (comm->isRoot())
//        UBLOG(logINFO, "  physicsEngineSolver->runTimestep() time = " << timer.stop() << " s");
//#endif
//}
//
//void DemCoProcessor::moveVfGeoObjects()
//{
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                    ->getSemiactive()) {
//                walberla::pe::RigidBody *peGeoObject = getPeGeoObject(
//                    std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                        ->getSystemID());
//                if (peGeoObject != nullptr) {
//                    std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                        ->setGeometry(peGeoObject);
//                    interactors[i]->moveGbObjectTo(physicsEngineGeometrieAdapters[i]->getPosition());
//                    std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                        ->setSemiactive(false);
//                } else {
//                    std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                        ->setInactive();
//                }
//            } else {
//                interactors[i]->moveGbObjectTo(physicsEngineGeometrieAdapters[i]->getPosition());
//            }
//        }
//    }
//}
//
//bool DemCoProcessor::isDemObjectInAABB(std::array<double, 6> AABB)
//{
//    bool result = false;
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            SPtr<GbObject3D> geoObject = interactors[i]->getGbObject3D();
//            std::array<double, 2> minMax1;
//            std::array<double, 2> minMax2;
//            std::array<double, 2> minMax3;
//            minMax1[0] = geoObject->getX1Minimum();
//            minMax2[0] = geoObject->getX2Minimum();
//            minMax3[0] = geoObject->getX3Minimum();
//            minMax1[1] = geoObject->getX1Maximum();
//            minMax2[1] = geoObject->getX2Maximum();
//            minMax3[1] = geoObject->getX3Maximum();
//
//            for (int x3 = 0; x3 < 2; x3++)
//                for (int x2 = 0; x2 < 2; x2++)
//                    for (int x1 = 0; x1 < 2; x1++) {
//                        result =
//                            result || (minMax1[x1] >= AABB[0] && minMax2[x2] >= AABB[1] && minMax3[x3] >= AABB[2] &&
//                                       minMax1[x1] <= AABB[3] && minMax2[x2] <= AABB[4] && minMax3[x3] <= AABB[5]);
//                    }
//        }
//    }
//
//    std::vector<int> values;
//    values.push_back((int)result);
//    std::vector<int> rvalues = comm->gather(values);
//
//    if (comm->isRoot()) {
//        for (int i = 0; i < (int)rvalues.size(); i++) {
//            result = result || (bool)rvalues[i];
//        }
//    }
//    int iresult = (int)result;
//    comm->broadcast(iresult);
//    result = (bool)iresult;
//
//    return result;
//}
//
//int DemCoProcessor::addSurfaceTriangleSet(std::vector<UbTupleFloat3> &nodes, std::vector<UbTupleInt3> &triangles)
//{
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            interactors[i]->getGbObject3D()->addSurfaceTriangleSet(nodes, triangles);
//        }
//    }
//    return (int)interactors.size();
//}
//
//void DemCoProcessor::getObjectsPropertiesVector(std::vector<double> &p)
//{
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            p.push_back(i);
//            p.push_back(interactors[i]->getGbObject3D()->getX1Centroid());
//            p.push_back(interactors[i]->getGbObject3D()->getX2Centroid());
//            p.push_back(interactors[i]->getGbObject3D()->getX3Centroid());
//            Vector3D v = physicsEngineGeometrieAdapters[i]->getLinearVelocity();
//            p.push_back(v[0]);
//            p.push_back(v[1]);
//            p.push_back(v[2]);
//        }
//    }
//}
//
//void DemCoProcessor::addPeGeo(walberla::pe::RigidBody *peGeo)
//{
//    auto geometry = getPeGeoAdapter(peGeo->getSystemID());
//    if (geometry != nullptr) {
//        geometry->setActive();
//        geometry->setGeometry(peGeo);
//        return;
//    } else
//        return;
//}
//
//void DemCoProcessor::removePeGeo(walberla::pe::RigidBody *peGeo)
//{
//    auto geometry = getPeGeoAdapter(peGeo->getSystemID());
//    if (geometry != nullptr) {
//        geometry->setSemiactive(true);
//    } else
//        throw UbException(UB_EXARGS, "PeGeo SystemId=" + UbSystem::toString(peGeo->getSystemID()) +
//                                         " is not matching geometry ID");
//}
//
//void DemCoProcessor::addPeShadowGeo(walberla::pe::RigidBody *peGeo)
//{
//    auto geometry = getPeGeoAdapter(peGeo->getSystemID());
//    if (geometry != nullptr) {
//        geometry->setActive();
//        geometry->setGeometry(peGeo);
//        return;
//    } else
//        throw UbException(UB_EXARGS,
//                          "PeGeo ID=" + UbSystem::toString(peGeo->getSystemID()) + " is not matching geometry ID");
//}
//
//void DemCoProcessor::removePeShadowGeo(walberla::pe::RigidBody *peGeo)
//{
//    auto geometry = getPeGeoAdapter(peGeo->getSystemID());
//
//    if (geometry != nullptr) {
//        geometry->setSemiactive(true);
//    } else
//        throw UbException(UB_EXARGS,
//                          "PeGeo ID=" + UbSystem::toString(peGeo->getSystemID()) + " is not matching geometry ID");
//}
//
//bool DemCoProcessor::isSpheresIntersection(double centerX1, double centerX2, double centerX3, double d)
//{
//    bool result = false;
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            SPtr<GbObject3D> sphere = interactors[i]->getGbObject3D();
//            result                  = result ||
//                     (sqrt(pow(sphere->getX1Centroid() - centerX1, 2.0) + pow(sphere->getX2Centroid() - centerX2, 2.0) +
//                           pow(sphere->getX3Centroid() - centerX3, 2.0)) <= d);
//        }
//    }
//    std::vector<int> values;
//    values.push_back((int)result);
//    std::vector<int> rvalues = comm->gather(values);
//
//    if (comm->isRoot()) {
//        for (int i = 0; i < (int)rvalues.size(); i++) {
//            result = result || (bool)rvalues[i];
//        }
//    }
//    int iresult = (int)result;
//    comm->broadcast(iresult);
//    result = (bool)iresult;
//
//    return result;
//}
//
//void DemCoProcessor::distributeIDs()
//{
//    std::vector<unsigned long long> peIDsSend;
//    std::vector<int> vfIDsSend;
//
//    for (int i = 0; i < interactors.size(); i++) {
//        if (std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])->isActive()) {
//            peIDsSend.push_back(
//                std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//                    ->getSystemID());
//            vfIDsSend.push_back(interactors[i]->getID());
//        }
//    }
//
//    std::vector<unsigned long long> peIDsRecv;
//    std::vector<int> vfIDsRecv;
//
//    comm->allGather(peIDsSend, peIDsRecv);
//    comm->allGather(vfIDsSend, vfIDsRecv);
//
//    std::map<int, unsigned long long> idMap;
//
//    for (int i = 0; i < peIDsRecv.size(); i++) {
//        idMap.insert(std::make_pair(vfIDsRecv[i], peIDsRecv[i]));
//    }
//
//    for (int i = 0; i < interactors.size(); i++) {
//        std::map<int, unsigned long long>::const_iterator it;
//        if ((it = idMap.find(interactors[i]->getID())) == idMap.end()) {
//            throw UbException(UB_EXARGS, "Interactor ID = " + UbSystem::toString(interactors[i]->getID()) +
//                                             " is invalid! The DEM object may be not in PE domain!");
//        }
//
//        std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])
//            ->setSystemID(it->second);
//
//        geoIdMap.insert(std::make_pair(
//            it->second, std::dynamic_pointer_cast<PePhysicsEngineGeometryAdapter>(physicsEngineGeometrieAdapters[i])));
//    }
//}
////////////////////////////////////////////////////////////////////////////
//void DemCoProcessor::setBlockVisitor(std::shared_ptr<BoundaryConditionsBlockVisitor> boundaryConditionsBlockVisitor)
//{
//    this->boundaryConditionsBlockVisitor = boundaryConditionsBlockVisitor;
//}
////////////////////////////////////////////////////////////////////////////
//walberla::pe::RigidBody *DemCoProcessor::getPeGeoObject(walberla::id_t id)
//{
//    std::shared_ptr<walberla::blockforest::BlockForest> forest =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getBlockForest();
//    std::shared_ptr<walberla::domain_decomposition::BlockDataID> storageId =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getStorageId();
//    std::shared_ptr<walberla::pe::BodyStorage> globalBodyStorage =
//        std::dynamic_pointer_cast<PePhysicsEngineSolverAdapter>(physicsEngineSolver)->getGlobalBodyStorage();
//
//    return walberla::pe::getBody(*globalBodyStorage, *forest, *storageId, id,
//                                 walberla::pe::StorageSelect::LOCAL | walberla::pe::StorageSelect::SHADOW);
//}
//////////////////////////////////////////////////////////////////////////////
//std::shared_ptr<PePhysicsEngineGeometryAdapter> DemCoProcessor::getPeGeoAdapter(unsigned long long systemId)
//{
//    std::map<unsigned long long, std::shared_ptr<PePhysicsEngineGeometryAdapter>>::const_iterator it;
//    if ((it = geoIdMap.find(systemId)) == geoIdMap.end()) {
//        return nullptr;
//    } else
//        return it->second;
//}

LiggghtsCoProcessor::LiggghtsCoProcessor()
{
}

LiggghtsCoProcessor::~LiggghtsCoProcessor()
{
}
