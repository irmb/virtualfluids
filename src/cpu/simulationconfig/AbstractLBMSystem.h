#ifndef VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H
#define VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H

#include <Interactors/Interactor3D.h>
#include <BoundaryConditions/BC.h>
#include <memory>

class AbstractLBMSystem {
public:
    virtual ~AbstractLBMSystem() = default;

    virtual int getNumberOfDirections() = 0;

    virtual std::shared_ptr<Interactor3D> makeInteractor() = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   int type) = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   std::shared_ptr<BC> bcAdapter, int type) = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   std::shared_ptr<BC> bcAdapter, int type, Interactor3D::Accuracy accuracy) = 0;

};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H
