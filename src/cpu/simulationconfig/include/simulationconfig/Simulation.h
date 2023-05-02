#ifndef VIRTUALFLUIDSBUILD_H
#define VIRTUALFLUIDSBUILD_H

#include <string>
#include <memory>
#include <set>

#include <mpi/Communicator.h>

#include <geometry3d/GbPoint3D.h>
#include <Interactors/Interactor3D.h>
#include <BoundaryConditions/BC.h>
#include <Visitors/BoundaryConditionsBlockVisitor.h>
#include <SimulationObservers/SimulationObserver.h>
#include <LBM/LBMUnitConverter.h>
#include "KernelFactory.h"
#include "AbstractLBMSystem.h"
#include "KernelConfigStructs.h"
#include "SimulationParameters.h"
#include "WriterConfiguration.h"



class CPUSimulation
{
private:
    KernelFactory kernelFactory = KernelFactory();

    std::shared_ptr<LBMKernel> lbmKernel;
    std::shared_ptr<AbstractLBMSystem> lbmSystem;
    std::shared_ptr<vf::mpi::Communicator> communicator;

    std::shared_ptr<Grid3D> grid;
    std::vector<std::shared_ptr<Interactor3D>> interactors;
    BoundaryConditionsBlockVisitor bcVisitor;
    std::set<std::shared_ptr<BC>> registeredAdapters;

    std::shared_ptr<LBMKernelConfiguration> kernelConfig;
    std::shared_ptr<RuntimeParameters> simulationParameters;
    std::shared_ptr<GridParameters> gridParameters;
    std::shared_ptr<PhysicalParameters> physicalParameters;

    WriterConfiguration &writerConfig = *(new WriterConfiguration());

public:
    explicit CPUSimulation();

    ~CPUSimulation();

    WriterConfiguration &getWriterConfig();

    void setWriterConfiguration(const WriterConfiguration &config);

    void setGridParameters(std::shared_ptr<GridParameters> parameters);

    void setPhysicalParameters(std::shared_ptr<PhysicalParameters> parameters);

    void setRuntimeParameters(std::shared_ptr<RuntimeParameters> parameters);

    void setKernelConfiguration(const std::shared_ptr<LBMKernelConfiguration> &kernel);

    void addObject(const std::shared_ptr<GbObject3D> &object, const std::shared_ptr<BC> &bcAdapter, int state,
                   const std::string &folderPath);

    void addBCAdapter(const std::shared_ptr<BC> &bcAdapter);

    void run();

private:
    bool isMainProcess();

    std::shared_ptr<GbObject3D> makeSimulationBoundingBox();

    void writeBlocksToFile() const;

    void writeBoundaryConditions() const;

    std::shared_ptr<SimulationObserver> makeMacroscopicQuantitiesCoProcessor(const std::shared_ptr<LBMUnitConverter> &converter,
                                                           const std::shared_ptr<UbScheduler> &visualizationScheduler) const;

    static std::shared_ptr<LBMUnitConverter> makeLBMUnitConverter();

    void setBlockSize(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;

    static void setBoundaryConditionProcessor(const std::shared_ptr<LBMKernel> &kernel);

    void generateBlockGrid(const std::shared_ptr<GbObject3D> &gridCube) const;

    void logSimulationData(const int &nodesInX1, const int &nodesInX2, const int &nodesInX3) const;


    void setKernelForcing(const std::shared_ptr<LBMKernel> &kernel, std::shared_ptr<LBMUnitConverter> &converter) const;

    void setConnectors();

    void initializeDistributions();

    std::shared_ptr<SimulationObserver> makeNupsCoProcessor() const;
};


#endif