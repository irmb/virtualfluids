//
// Created by Sven Marcus on 20.08.20.
//

#ifndef VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H
#define VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H

#include <Interactors/Interactor3D.h>
#include <BoundaryConditions/BCAdapter.h>
#include <memory>

class AbstractLBMSystem {
public:
    virtual int getNumberOfDirections() = 0;

    virtual std::shared_ptr<Interactor3D> makeInteractor() = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   int type) = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   std::shared_ptr<BCAdapter> bcAdapter, int type) = 0;

    virtual std::shared_ptr<Interactor3D>
    makeInteractor(std::shared_ptr<GbObject3D> object, std::shared_ptr<Grid3D> grid,
                   std::shared_ptr<BCAdapter> bcAdapter, int type, Interactor3D::Accuracy accuracy) = 0;

};


#endif //VIRTUALFLUIDSPYTHONBINDINGS_ABSTRACTLBMSYSTEM_H
