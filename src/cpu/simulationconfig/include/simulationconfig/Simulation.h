#ifndef VIRTUALFLUIDSBUILD_H
#define VIRTUALFLUIDSBUILD_H

#include <string>
#include <memory>
#include <set>
#include <geometry3d/GbPoint3D.h>
#include <Interactors/Interactor3D.h>
#include <BoundaryConditions/BCAdapter.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <CoProcessors/CoProcessor.h>
#include <LBM/LBMUnitConverter.h>
#include <Parallel/Communicator.h>
#include "KernelFactory.h"
#include "AbstractLBMSystem.h"
#include "KernelConfigStructs.h"
#include "SimulationParameters.h"
#include "WriterConfiguration.h"


class Simulation {
private:
    KernelFactory kernelFactory = KernelFactory();

    std::shared_ptr<LBMKernel> lbmKernel;
    std::shared_ptr<AbstractLBMSystem> lbmSystem;
    std::shared_ptr<Communicator> communicator;

    std::shared_ptr<Grid3D> grid;
    std::vector<std::shared_ptr<Interactor3D>> interactors;
    BoundaryConditionsBlockVisitor bcVisitor;
    std::set<std::shared_ptr<BCAdapter>> registeredAdapters;

    std::shared_ptr<LBMKernelConfiguration> kernelConfig;
    std::shared_ptr<RuntimeParameters> simulationParameters;
    std::shared_ptr<GridParameters> gridParameters;
    std::shared_ptr<PhysicalParameters> physicalParameters;

    WriterConfiguration &writerConfig = *(new WriterConfiguration());

public:
    explicit Simulation();

    ~Simulation();

    WriterConfiguration &getWriterConfig();

    void setWriterConfiguration(const WriterConfiguration &config);

    void setGridParameters(std::shared_ptr<GridParameters> parameters);

    void setPhysicalParameters(std::shared_ptr<PhysicalParameters> parameters);

    void setRuntimeParameters(std::shared_ptr<RuntimeParameters> parameters);

    void setKernelConfiguration(const std::shared_ptr<LBMKernelConfiguration> &kernel);

    void addObject(const std::shared_ptr<GbObject3D> &object, const std::shared_ptr<BCAdapter> &bcAdapter, int state,
                   const std::string &folderPath);

    void addBCAdapter(const std::shared_ptr<BCAdapter> &bcAdapter);

    void run();

private:
    std::shared_ptr<GbObject3D> makeSimulationBoundingBox(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    void writeBlocks() const;

    void writeBoundaryConditions() const;

    std::shared_ptr<CoProcessor> makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                           const std::shared_ptr<UbScheduler> &visSch) const;

    static std::shared_ptr<LBMUnitConverter> makeLBMUnitConverter();

    void setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    static void setBoundaryConditionProcessor(const std::shared_ptr<LBMKernel> &kernel);

    void generateBlockGrid(const std::shared_ptr<GbObject3D> &gridCube) const;

    void logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;


    void setKernelForcing(const std::shared_ptr<LBMKernel> &kernel, std::shared_ptr<LBMUnitConverter> &converter) const;
};

#endif