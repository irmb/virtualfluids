#ifndef VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H
#define VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H


#include "AbstractLBMSystem.h"

class D3Q27LBMSystem : public AbstractLBMSystem {
public:
    D3Q27LBMSystem() = default;

    int getNumberOfDirections() override;

    std::shared_ptr<Interactor3D> makeInteractor() override;

    std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid, int type) override;

    std::shared_ptr<Interactor3D> makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                                                 std::shared_ptr<BC> bcAdapter, int type) override;

    std::shared_ptr<Interactor3D> makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                                                 std::shared_ptr<BC> bcAdapter, int type,
                                                 Interactor3D::Accuracy accuracy) override;
};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_D3Q27LBMSYSTEM_H
