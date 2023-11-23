#ifdef VF_CATALYST

#include "InSituCatalystSimulationObserver.h"
#include <D3Q27ETBCSet.h>
#include <LBMKernel.h>
#include <string>
#include <vector>

#include <vtkCellType.h>
#include <vtkPointData.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <fstream>
#include <iostream>

using namespace std;

InSituCatalystSimulationObserver::InSituCatalystSimulationObserver() {}
//////////////////////////////////////////////////////////////////////////
InSituCatalystSimulationObserver::InSituCatalystSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, std::string script)
    : SimulationObserver(grid, s)
{
    gridRank     = vf::parallel::Communicator::getInstance()->getProcessID();
    minInitLevel = this->grid->getCoarsestInitializedLevel();
    maxInitLevel = this->grid->getFinestInitializedLevel();

    blockVector.resize(maxInitLevel + 1);

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        grid->getBlocks(level, gridRank, true, blockVector[level]);
    }

    Processor = vtkSmartPointer<vtkCPProcessor>::New();
    Processor->Initialize();

    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize(script.c_str());
    Processor->AddPipeline(pipeline.GetPointer());

    buildVTKGrid();
}
//////////////////////////////////////////////////////////////////////////
InSituCatalystSimulationObserver::~InSituCatalystSimulationObserver() {}
//////////////////////////////////////////////////////////////////////////
void InSituCatalystSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);

    UBLOG(logDEBUG3, "InSituCatalystSimulationObserver::update:" << step);
}
//////////////////////////////////////////////////////////////////////////
void InSituCatalystSimulationObserver::collectData(real step)
{
    unsigned int istep = static_cast<int>(step);

    vtkNew<vtkCPDataDescription> dataDescription;
    dataDescription->AddInput("input");
    dataDescription->SetTimeData(step, istep);

    if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0) {

        index = 0;

        for (int level = minInitLevel; level <= maxInitLevel; level++) {
            for (SPtr<Block3D> block : blockVector[level]) {
                if (block) {
                    addData(block);
                }
            }
        }

        vtkDoubleArray *rho = vtkDoubleArray::SafeDownCast(unstructuredGrid->GetPointData()->GetArray("Rho"));
        rho->SetArray(&rhoArray[0], static_cast<vtkIdType>(rhoArray.size()), 1);

        vtkDoubleArray *vx1 = vtkDoubleArray::SafeDownCast(unstructuredGrid->GetPointData()->GetArray("Vx"));
        vx1->SetArray(&vx1Array[0], static_cast<vtkIdType>(vx1Array.size()), 1);

        vtkDoubleArray *vx2 = vtkDoubleArray::SafeDownCast(unstructuredGrid->GetPointData()->GetArray("Vy"));
        vx2->SetArray(&vx2Array[0], static_cast<vtkIdType>(vx2Array.size()), 1);

        vtkDoubleArray *vx3 = vtkDoubleArray::SafeDownCast(unstructuredGrid->GetPointData()->GetArray("Vz"));
        vx3->SetArray(&vx3Array[0], static_cast<vtkIdType>(vx3Array.size()), 1);

        dataDescription->GetInputDescriptionByName("input")->SetGrid(unstructuredGrid);
        Processor->CoProcess(dataDescription.GetPointer());
    }

    UBLOG(logINFO, "InSituCatalystSimulationObserver step: " << istep);
}
//////////////////////////////////////////////////////////////////////////
void InSituCatalystSimulationObserver::addData(SPtr<Block3D> block)
{
    UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
    real dx                   = grid->getDeltaX(block);

    SPtr<LBMKernel> kernel                  = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)(distributions->getNX1());
    int maxX2 = (int)(distributions->getNX2());
    int maxX3 = (int)(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

    for (size_t ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (size_t ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (size_t ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    distributions->getDistribution(f, ix1, ix2, ix3);
                    calcMacros(f, rho, vx1, vx2, vx3);
                    real press = D3Q27System::calcPress(f, rho, vx1, vx2, vx3);

                    if (UbMath::isNaN(rho) || UbMath::isInfinity(rho))
                        UB_THROW(UbException(
                            UB_EXARGS, "rho is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // rho=999.0;
                    if (UbMath::isNaN(press) || UbMath::isInfinity(press))
                        UB_THROW(UbException(
                            UB_EXARGS, "press is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // press=999.0;
                    if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx1=999.0;
                    if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx2=999.0;
                    if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx3=999.0;

                    rhoArray[index] = rho;
                    vx1Array[index] = vx1;
                    vx2Array[index] = vx2;
                    vx3Array[index] = vx3;
                    index++;
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////
void InSituCatalystSimulationObserver::buildVTKGrid()
{
    unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    points           = vtkPoints::New();
    unstructuredGrid->SetPoints(points);
    arrays[0] = vtkSmartPointer<vtkDoubleArray>::New();
    arrays[0]->SetNumberOfComponents(1);
    arrays[0]->SetName("Rho");
    arrays[1] = vtkSmartPointer<vtkDoubleArray>::New();
    arrays[1]->SetNumberOfComponents(1);
    arrays[1]->SetName("Vx");
    arrays[2] = vtkSmartPointer<vtkDoubleArray>::New();
    arrays[2]->SetNumberOfComponents(1);
    arrays[2]->SetName("Vy");
    arrays[3] = vtkSmartPointer<vtkDoubleArray>::New();
    arrays[3]->SetNumberOfComponents(1);
    arrays[3]->SetName("Vz");

    numOfPoints = 0;

    for (int level = minInitLevel; level <= maxInitLevel; level++) {
        for (SPtr<Block3D> block : blockVector[level]) {
            if (block) {
                addVTKGridData(block);
            }
        }
    }

    unstructuredGrid->GetPointData()->AddArray(arrays[0]);
    unstructuredGrid->GetPointData()->AddArray(arrays[1]);
    unstructuredGrid->GetPointData()->AddArray(arrays[2]);
    unstructuredGrid->GetPointData()->AddArray(arrays[3]);
    unstructuredGrid->GetPointData()->SetScalars(arrays[1]);

    rhoArray.resize(numOfPoints);
    vx1Array.resize(numOfPoints);
    vx2Array.resize(numOfPoints);
    vx3Array.resize(numOfPoints);
}
//////////////////////////////////////////////////////////////////////////
void InSituCatalystSimulationObserver::addVTKGridData(SPtr<Block3D> block)
{
    UbTupleDouble3 org          = grid->getBlockWorldCoordinates(block);
    UbTupleDouble3 blockLengths = grid->getBlockLengths(block);
    UbTupleDouble3 nodeOffset   = grid->getNodeOffset(block);
    real dx                   = grid->getDeltaX(block);

    SPtr<LBMKernel> kernel                  = block->getKernel();
    SPtr<BCArray3D> bcArray                 = kernel->getBCSet()->getBCArray();
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();
    real f[D3Q27System::ENDF + 1];
    real vx1, vx2, vx3, rho;

    // knotennummerierung faengt immer bei 0 an!
    int SWB, SEB, NEB, NWB, SWT, SET, NET, NWT;

    // Funktionszeiger
    // typedef void(*CalcMacrosFct)(const real* const& /*feq[27]*/, real& /*(d)rho*/, real& /*vx1*/, real&
    // /*vx2*/, real& /*vx3*/);

    // CalcMacrosFct calcMacros = NULL;

    if (block->getKernel()->getCompressible()) {
        calcMacros = &D3Q27System::calcCompMacroscopicValues;
    } else {
        calcMacros = &D3Q27System::calcIncompMacroscopicValues;
    }

    int minX1 = 0;
    int minX2 = 0;
    int minX3 = 0;

    int maxX1 = (int)(distributions->getNX1());
    int maxX2 = (int)(distributions->getNX2());
    int maxX3 = (int)(distributions->getNX3());

    // nummern vergeben und node vector erstellen + daten sammeln
    CbArray3D<int> nodeNumbers((int)maxX1, (int)maxX2, (int)maxX3, -1);
    maxX1 -= 2;
    maxX2 -= 2;
    maxX3 -= 2;

    SPtr<BoundaryConditions> bcPtr;
    int nr = points->GetNumberOfPoints();

    real x[3];

    for (size_t ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (size_t ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (size_t ix1 = minX1; ix1 <= maxX1; ix1++) {
                if (!bcArray->isUndefined(ix1, ix2, ix3) && !bcArray->isSolid(ix1, ix2, ix3)) {
                    x[0] = real(val<1>(org) - val<1>(nodeOffset) + ix1 * dx);
                    x[1] = real(val<2>(org) - val<2>(nodeOffset) + ix2 * dx);
                    x[2] = real(val<3>(org) - val<3>(nodeOffset) + ix3 * dx);

                    points->InsertPoint((vtkIdType)nr, x);

                    nodeNumbers(ix1, ix2, ix3) = nr++;

                    distributions->getDistribution(f, ix1, ix2, ix3);
                    calcMacros(f, rho, vx1, vx2, vx3);
                    real press = D3Q27System::calcPress(f, rho, vx1, vx2, vx3);

                    if (UbMath::isNaN(rho) || UbMath::isInfinity(rho))
                        UB_THROW(UbException(
                            UB_EXARGS, "rho is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // rho=999.0;
                    if (UbMath::isNaN(press) || UbMath::isInfinity(press))
                        UB_THROW(UbException(
                            UB_EXARGS, "press is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // press=999.0;
                    if (UbMath::isNaN(vx1) || UbMath::isInfinity(vx1))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx1 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx1=999.0;
                    if (UbMath::isNaN(vx2) || UbMath::isInfinity(vx2))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx2 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx2=999.0;
                    if (UbMath::isNaN(vx3) || UbMath::isInfinity(vx3))
                        UB_THROW(UbException(
                            UB_EXARGS, "vx3 is not a number (nan or -1.#IND) or infinity number -1.#INF in block=" +
                                           block->toString() + ", node=" + UbSystem::toString(ix1) + "," +
                                           UbSystem::toString(ix2) + "," + UbSystem::toString(ix3)));
                    // vx3=999.0;

                    arrays[0]->InsertNextValue(rho);
                    arrays[1]->InsertNextValue(vx1);
                    arrays[2]->InsertNextValue(vx2);
                    arrays[3]->InsertNextValue(vx3);
                    numOfPoints++;
                }
            }
        }
    }
    maxX1 -= 1;
    maxX2 -= 1;
    maxX3 -= 1;

    vtkIdType ptIds[8];
    // cell vector erstellen
    for (int ix3 = minX3; ix3 <= maxX3; ix3++) {
        for (int ix2 = minX2; ix2 <= maxX2; ix2++) {
            for (int ix1 = minX1; ix1 <= maxX1; ix1++) {
                if ((ptIds[0] = SWB = nodeNumbers(ix1, ix2, ix3)) >= 0 &&
                    (ptIds[1] = SEB = nodeNumbers(ix1 + 1, ix2, ix3)) >= 0 &&
                    (ptIds[2] = NWB = nodeNumbers(ix1, ix2 + 1, ix3)) >= 0 &&
                    (ptIds[3] = NEB = nodeNumbers(ix1 + 1, ix2 + 1, ix3)) >= 0 &&
                    (ptIds[4] = SWT = nodeNumbers(ix1, ix2, ix3 + 1)) >= 0 &&
                    (ptIds[5] = SET = nodeNumbers(ix1 + 1, ix2, ix3 + 1)) >= 0 &&
                    (ptIds[6] = NWT = nodeNumbers(ix1, ix2 + 1, ix3 + 1)) >= 0 &&
                    (ptIds[7] = NET = nodeNumbers(ix1 + 1, ix2 + 1, ix3 + 1)) >= 0) {
                    unstructuredGrid->InsertNextCell((int)VTK_VOXEL, (vtkIdType)8, ptIds);
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////

#endif
