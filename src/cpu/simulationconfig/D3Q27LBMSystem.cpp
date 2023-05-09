#include "simulationconfig/D3Q27LBMSystem.h"
#include <Interactors/D3Q27Interactor.h>
#include <LBM/D3Q27System.h>

int D3Q27LBMSystem::getNumberOfDirections()
{
    return D3Q27System::ENDDIR;
}

std::shared_ptr<Interactor3D> D3Q27LBMSystem::makeInteractor()
{
    return std::shared_ptr<Interactor3D>(new D3Q27Interactor());
}

std::shared_ptr<Interactor3D>
D3Q27LBMSystem::makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid, int type)
{
    return std::shared_ptr<Interactor3D>(new D3Q27Interactor(object, grid, type));
}

std::shared_ptr<Interactor3D>
D3Q27LBMSystem::makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                               std::shared_ptr<BC> bcAdapter, int type)
{
    return std::shared_ptr<Interactor3D>(new D3Q27Interactor(object, grid, bcAdapter, type));
}

std::shared_ptr<Interactor3D>
D3Q27LBMSystem::makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                               std::shared_ptr<BC> bcAdapter, int type, Interactor3D::Accuracy accuracy)
{
    return std::shared_ptr<Interactor3D>(new D3Q27Interactor(object, grid, bcAdapter, type, accuracy));
}
