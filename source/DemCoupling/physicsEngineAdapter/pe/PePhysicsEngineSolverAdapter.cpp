#include "PePhysicsEngineSolverAdapter.h"

#include <exception>

#include <pe/basic.h>
#include "PeAdapter.h"
#include "PePhysicsEngineGeometryAdapter.h"
#include "PePhysicsEngineMaterialAdapter.h"
#include "PeLoadBalancerAdapter.h"
#include "Communicator.h"
#include "UbLogger.h"
#include <boost/tuple/tuple.hpp>
#include "UbException.h"
#include "UbSystem.h"
#include <memory>

typedef boost::tuple<walberla::pe::Sphere, walberla::pe::Plane> BodyTypeTuple;


PePhysicsEngineSolverAdapter::PePhysicsEngineSolverAdapter(std::shared_ptr<PeParameter> peParameter, std::shared_ptr<PeLoadBalancerAdapter> loadBalancer) : peParameter(peParameter), loadBalancer(loadBalancer)
{
    this->initalizePeEnvironment();
}

void PePhysicsEngineSolverAdapter::initalizePeEnvironment()
{
    this->initialPeBodyStorage();
    this->initialPeBlockForest();
    this->initalBlockData();
    this->initalPeIntegrator();
    this->executePeBodyTypeTuple();
    this->initialPeChannel();
}


std::shared_ptr<PhysicsEngineGeometryAdapter> PePhysicsEngineSolverAdapter::createPhysicsEngineGeometryAdapter(int id, const Vector3D& position, double radius, std::shared_ptr<PhysicsEngineMaterialAdapter> material) const
{
    const std::shared_ptr<PePhysicsEngineMaterialAdapter> peMaterial = std::dynamic_pointer_cast<PePhysicsEngineMaterialAdapter>(material);
    std::shared_ptr<PePhysicsEngineGeometryAdapter> peGeometryAdapter(new PePhysicsEngineGeometryAdapter());

    //UBLOG(logINFO, "PePhysicsEngineSolverAdapter::createSphere():start");
    walberla::pe::GeomID peGeometry = createSphere(*globalBodyStorage, *forest, *storageId, id, PeConverter::convert(position), radius, peMaterial->getPeMaterial());
    //UBLOG(logINFO, "PePhysicsEngineSolverAdapter::createSphere():end");

    if (peGeometry)
    {
       peGeometryAdapter->setId(id);
       peGeometryAdapter->setSystemID(peGeometry->getSystemID());
       //unsigned long long sid = peGeometryAdapter->getSystemID();
       //geoIdMap.insert(std::make_pair(sid, peGeometryAdapter));
       //peGeometryAdapter->counter++;
       
       //if(peGeometry->getSystemID() != id+1)
       //   UB_THROW(UbException(UB_EXARGS, "id="+UbSystem::toString(id)+" does not match pe::SystemId="+UbSystem::toString(peGeometry->getSystemID())));

       peGeometryAdapter->setActive();
       peGeometryAdapter->setGeometry(peGeometry);
       return peGeometryAdapter;
    }
    else
    {
       peGeometryAdapter->setId(id);
       peGeometryAdapter->setInactive();
       return peGeometryAdapter;
    }

    walberla::pe::syncNextNeighbors<BodyTypeTuple>(*forest, *storageId);
    //walberla::pe::syncShadowOwners<BodyTypeTuple>(*forest, *storageId);
}

void PePhysicsEngineSolverAdapter::runTimestep(double step)
{
    cr->timestep(walberla::real_c(step));
    walberla::pe::syncNextNeighbors<BodyTypeTuple>(*forest, *storageId);
    //walberla::pe::syncShadowOwners<BodyTypeTuple>(*forest, *storageId);
}



void PePhysicsEngineSolverAdapter::initialPeBodyStorage()
{
    globalBodyStorage = std::make_shared<walberla::pe::BodyStorage>();
}

void PePhysicsEngineSolverAdapter::initialPeBlockForest()
{

   //walberla::SetupBlockForest sforest = walberla::blockforest::createUniformBlockGrid(walberla::AABB(peParameter->simulationDomain[0], peParameter->simulationDomain[1], peParameter->simulationDomain[2],
   //   peParameter->simulationDomain[3], peParameter->simulationDomain[4], peParameter->simulationDomain[5]), // simulationDomain
   //   walberla::uint_t(val<1>(peParameter->numberOfBlocks)), walberla::uint_t(val<2>(peParameter->numberOfBlocks)), walberla::uint_t(val<3>(peParameter->numberOfBlocks)),walberla::uint_t(10),walberla::uint_t(10),walberla::uint_t(10), 5.0,false);
    walberla::SetupBlockForest sforest;
    //sforest.addWorkloadMemorySUIDAssignmentFunction( uniformWorkloadAndMemoryAssignment );
    sforest.init(walberla::AABB(peParameter->simulationDomain[0], peParameter->simulationDomain[1], peParameter->simulationDomain[2],
       peParameter->simulationDomain[3], peParameter->simulationDomain[4], peParameter->simulationDomain[5]), // simulationDomain
       walberla::uint_t(val<1>(peParameter->numberOfBlocks)), walberla::uint_t(val<2>(peParameter->numberOfBlocks)), walberla::uint_t(val<3>(peParameter->numberOfBlocks)), // blocks in each direction
       val<1>(peParameter->isPeriodic), val<2>(peParameter->isPeriodic), val<3>(peParameter->isPeriodic));
    sforest.balanceLoad(*loadBalancer.get(), loadBalancer->getNumberOfProcesses());
    forest = std::shared_ptr< walberla::blockforest::BlockForest >( new walberla::blockforest::BlockForest( walberla::uint_c( loadBalancer->getRank() ), sforest) );

     auto mpiManager = walberla::MPIManager::instance();
     mpiManager->useWorldComm();
    if (!forest)
       throw std::runtime_error("No PE BlockForest created ... ");
}

void PePhysicsEngineSolverAdapter::initalBlockData()
{
    storageId = std::make_shared<walberla::domain_decomposition::BlockDataID>
    (
        forest->addBlockData(walberla::pe::createStorageDataHandling<BodyTypeTuple>(), "Storage")
    );
}

void PePhysicsEngineSolverAdapter::initalPeIntegrator()
{
    auto ccdID = forest->addBlockData(walberla::pe::ccd::createHashGridsDataHandling(globalBodyStorage, *storageId), "CCD");
    auto fcdID = forest->addBlockData(walberla::pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, walberla::pe::fcd::AnalyticCollideFunctor>(), "FCD");

    cr = std::make_shared<walberla::pe::cr::HardContactSemiImplicitTimesteppingSolvers>(globalBodyStorage, forest, *storageId, ccdID, fcdID);
    cr->setMaxIterations(peParameter->maxPeIterations);
    cr->setRelaxationModel(walberla::pe::cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling);
    cr->setRelaxationParameter(walberla::real_t(peParameter->relaxationParameter));
    cr->setGlobalLinearAcceleration(PeConverter::convert(peParameter->globalLinearAcceleration));
}

void PePhysicsEngineSolverAdapter::executePeBodyTypeTuple()
{
    walberla::pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
}

void PePhysicsEngineSolverAdapter::initialPeChannel() const
{
    const walberla::pe::MaterialID material = peParameter->planes->getPeMaterial();

    auto simulationDomain = forest->getDomain();

    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(1, 0, 0), simulationDomain.minCorner(), material);
    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(-1, 0, 0), simulationDomain.maxCorner(), material);
    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 1, 0), simulationDomain.minCorner(), material);
    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, -1, 0), simulationDomain.maxCorner(), material);
    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 0, 1), simulationDomain.minCorner(), material);
    //createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 0, -1), simulationDomain.maxCorner(), material);

    Vector3D minOffset = peParameter->minOffset;
    Vector3D maxOffset = peParameter->maxOffset;

    walberla::pe::Vec3 minX1_Offset( minOffset.X1(), 0, 0);
    walberla::pe::Vec3 maxX1_Offset( maxOffset.X1(), 0, 0);
    walberla::pe::Vec3 minX2_Offset( 0, minOffset.X2(), 0);
    walberla::pe::Vec3 maxX2_Offset( 0, maxOffset.X2(), 0);
    walberla::pe::Vec3 minX3_Offset( 0, 0, minOffset.X3());
    walberla::pe::Vec3 maxX3_Offset( 0, 0, maxOffset.X3());

    walberla::pe::Vec3 minCorner = simulationDomain.minCorner();
    walberla::pe::Vec3 maxCorner = simulationDomain.maxCorner();

    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3( 1, 0, 0), minCorner + minX1_Offset, material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(-1, 0, 0), maxCorner + maxX1_Offset, material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3( 0, 1, 0), minCorner + minX2_Offset, material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3( 0,-1, 0), maxCorner + maxX2_Offset, material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3( 0, 0, 1), minCorner + minX3_Offset, material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3( 0, 0,-1), maxCorner + maxX3_Offset, material);
}

std::shared_ptr< walberla::blockforest::BlockForest > PePhysicsEngineSolverAdapter::getForest()
{
   return forest;
}

void PePhysicsEngineSolverAdapter::saveToFile(const std::string & path)
{
   //forest->saveToFile(path+"SerializeDeserialize.sbf");
   //forest->saveBlockData("SerializeDeserialize.dump", *storageId.get());
}

void PePhysicsEngineSolverAdapter::loadFromFile(const std::string & path)
{
   //forest = std::make_shared< walberla::blockforest::BlockForest >( walberla::uint_c( walberla::MPIManager::instance()->rank() ), path+"SerializeDeserialize.sbf", true, false );
   //storageId = std::make_shared< walberla::domain_decomposition::BlockDataID >(forest->loadBlockData(path+"SerializeDeserialize.dump", walberla::pe::createStorageDataHandling<BodyTypeTuple>(), "Storage"));
   //
   //auto ccdID = forest->addBlockData(walberla::pe::ccd::createHashGridsDataHandling(globalBodyStorage, *storageId), "CCD");
   //auto fcdID = forest->addBlockData(walberla::pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, walberla::pe::fcd::AnalyticCollideFunctor>(), "FCD");

   //cr = std::make_shared<walberla::pe::cr::HardContactSemiImplicitTimesteppingSolvers>(globalBodyStorage, forest, *storageId, ccdID, fcdID);
   //cr->setMaxIterations(peParameter->maxPeIterations);
   //cr->setRelaxationModel(walberla::pe::cr::HardContactSemiImplicitTimesteppingSolvers::ApproximateInelasticCoulombContactByDecoupling);
   //cr->setRelaxationParameter(walberla::real_t(peParameter->relaxationParameter));
   //cr->setGlobalLinearAcceleration(PeConverter::convert(peParameter->globalLinearAcceleration));

   //for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
   //{
   //   walberla::pe::ccd::ICCD* ccd = blockIt->getData< walberla::pe::ccd::ICCD >(ccdID);
   //   ccd->reloadBodies();
   //}
}

std::shared_ptr<walberla::blockforest::BlockForest> PePhysicsEngineSolverAdapter::getBlockForest()
{
   return forest;
}

std::shared_ptr<walberla::domain_decomposition::BlockDataID> PePhysicsEngineSolverAdapter::getStorageId()
{
   return storageId;
}

std::shared_ptr<walberla::pe::BodyStorage> PePhysicsEngineSolverAdapter::getGlobalBodyStorage()
{
   return globalBodyStorage;
}
//////////////////////////////////////////////////////////////////////////
//std::shared_ptr<PePhysicsEngineGeometryAdapter> PePhysicsEngineSolverAdapter::getPeGeoAdapter(unsigned long long systemId)
//{
//   std::map< unsigned long long, std::shared_ptr<PePhysicsEngineGeometryAdapter> >::const_iterator it;
//   if ((it=geoIdMap.find(systemId)) == geoIdMap.end())
//   {
//      //throw UbException(UB_EXARGS, "pe SystemId can not be found!");
//      return nullptr;
//   }
//   else
//      return it->second;
//}