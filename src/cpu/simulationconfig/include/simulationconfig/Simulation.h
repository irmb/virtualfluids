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
#include "WriterConfig.h"


class Simulation {
private:
    KernelFactory kernelFactory = KernelFactory();

    std::shared_ptr<LBMKernel> lbmKernel;
    std::shared_ptr<AbstractLBMSystem> lbmSystem;
    std::shared_ptr<Communicator> communicator;

    SPtr<Grid3D> grid;
    std::vector<SPtr<Interactor3D>> interactors;
    BoundaryConditionsBlockVisitor bcVisitor;
    std::set<SPtr<BCAdapter>> registeredAdapters;

    SPtr<LBMKernelConfig> kernelConfig;
    SPtr<SimulationParameters> simulationParameters;
    SPtr<GridParameters> gridParameters;
    SPtr<PhysicalParameters> physicalParameters;

    WriterConfig &writerConfig = *(new WriterConfig());

public:
    explicit Simulation();

    ~Simulation();

    WriterConfig &getWriterConfig();

    void setWriterConfig(const WriterConfig &config);

    void setGridParameters(SPtr<GridParameters> parameters);

    void setPhysicalParameters(SPtr<PhysicalParameters> parameters);

    void setSimulationParameters(SPtr<SimulationParameters> parameters);

    void setKernelConfig(const SPtr<LBMKernelConfig> &kernel);

    void addObject(const SPtr<GbObject3D> &object, const SPtr<BCAdapter> &bcAdapter, int state,
                   const std::string &folderPath);

    void addBCAdapter(const SPtr<BCAdapter> &bcAdapter);

    void run();

private:
    SPtr<GbObject3D> makeSimulationBoundingBox(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    void writeBlocks() const;

    void writeBoundaryConditions() const;

    SPtr<CoProcessor> makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                           const SPtr<UbScheduler> &visSch) const;

    static std::shared_ptr<LBMUnitConverter> makeLBMUnitConverter();

    void setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    static void setBoundaryConditionProcessor(const SPtr<LBMKernel> &kernel);

    void generateBlockGrid(const SPtr<GbObject3D> &gridCube) const;

    void logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;


    void setKernelForcing(const SPtr<LBMKernel> &kernel, std::shared_ptr<LBMUnitConverter> &converter) const;
};

#endif