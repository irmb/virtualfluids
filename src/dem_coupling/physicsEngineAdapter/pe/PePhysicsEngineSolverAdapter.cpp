#include "PePhysicsEngineSolverAdapter.h"

#include <exception>

#include <pe/basic.h>
#include "PeAdapter.h"
#include "PePhysicsEngineGeometryAdapter.h"
#include "PePhysicsEngineMaterialAdapter.h"

typedef boost::tuple<walberla::pe::Sphere, walberla::pe::Plane> BodyTypeTuple;


PePhysicsEngineSolverAdapter::PePhysicsEngineSolverAdapter(std::shared_ptr<PeParameter> peParameter) : peParameter(peParameter)
{
    this->initalizePeEnvironment();
}

void PePhysicsEngineSolverAdapter::initalizePeEnvironment()
{
    this->initialPeBodyStorage();
    this->initalPeBlockForest();
    this->initalBlockData();
    this->initalPeIntegrator();
    this->executePeBodyTypeTuple();
    this->initalPeChannel();
}


std::shared_ptr<PhysicsEngineGeometryAdapter> PePhysicsEngineSolverAdapter::createPhysicsEngineGeometryAdapter(int id, const Vector3D& position, double radius, std::shared_ptr<PhysicsEngineMaterialAdapter> material) const
{
    const std::shared_ptr<PePhysicsEngineMaterialAdapter> peMaterial = std::dynamic_pointer_cast<PePhysicsEngineMaterialAdapter>(material);

    walberla::pe::GeomID peGeometry = createSphere(*globalBodyStorage, *forest, *storageId, id, PeConverter::convert(position), radius, peMaterial->getPeMaterial());

    return std::static_pointer_cast<PhysicsEngineGeometryAdapter>(
        std::make_shared<PePhysicsEngineGeometryAdapter>(peGeometry));
}

void PePhysicsEngineSolverAdapter::runTimestep(double step)
{
    cr->timestep(walberla::real_c(step));
    walberla::pe::syncNextNeighbors<BodyTypeTuple>(*forest, *storageId);
}



void PePhysicsEngineSolverAdapter::initialPeBodyStorage()
{
    globalBodyStorage = boost::make_shared<walberla::pe::BodyStorage>();
}

void PePhysicsEngineSolverAdapter::initalPeBlockForest()
{
    forest = walberla::pe::createBlockForest
    (
        walberla::AABB(0, 0, 0, val<1>(peParameter->simulationDomain), val<2>(peParameter->simulationDomain), val<3>(peParameter->simulationDomain)), // simulationDomain
        walberla::Vector3<walberla::uint_t>(val<1>(peParameter->numberOfBlocks), val<2>(peParameter->numberOfBlocks), val<3>(peParameter->numberOfBlocks)), // blocks in each direction
        walberla::Vector3<bool>(val<1>(peParameter->isPeriodic), val<2>(peParameter->isPeriodic), val<3>(peParameter->isPeriodic)) // periodicity
    ); 

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

void PePhysicsEngineSolverAdapter::initalPeChannel() const
{
    const walberla::pe::MaterialID material = peParameter->planes->getPeMaterial();

    auto simulationDomain = forest->getDomain();
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(1, 0, 0), simulationDomain.minCorner(), material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(-1, 0, 0), simulationDomain.maxCorner(), material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 1, 0), simulationDomain.minCorner(), material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, -1, 0), simulationDomain.maxCorner(), material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 0, 1), simulationDomain.minCorner(), material);
    createPlane(*globalBodyStorage, 0, walberla::pe::Vec3(0, 0, -1), simulationDomain.maxCorner(), material);
}

walberla::pe::RigidBody* PePhysicsEngineSolverAdapter::getPeGeoObject(int id)
{
    for (auto blockIt = forest->begin(); blockIt != forest->end(); ++blockIt)
    {
        for (auto bodyIt = walberla::pe::BodyIterator::begin(*blockIt, *storageId.get()); bodyIt != walberla::pe::BodyIterator::end(); ++bodyIt)
        {
            if(bodyIt->getID() == id)
            {
                walberla::pe::RigidBody* geo = *(bodyIt);
                return geo;
            }
        }
    }
}

void PePhysicsEngineSolverAdapter::updateGeometry(std::shared_ptr<PhysicsEngineGeometryAdapter> geometryAdapter)
{
    auto geometry = std::static_pointer_cast<PePhysicsEngineGeometryAdapter>(geometryAdapter);
    const int id = geometry->getId();
    geometry->setGeometry(this->getPeGeoObject(id));
}
