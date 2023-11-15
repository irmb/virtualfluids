#ifndef GPU_K17CompressibleNavierStokes_H
#define GPU_K17CompressibleNavierStokes_H

#include "Kernel/KernelImp.h"

#include <lbm/collision/TurbulentViscosity.h>

class Parameter;

template <vf::lbm::TurbulenceModel turbulenceModel>
class K17CompressibleNavierStokes : public KernelImp
{
public:
    static std::shared_ptr<K17CompressibleNavierStokes<turbulenceModel>> getNewInstance(std::shared_ptr<Parameter> para,
                                                                                        int level);
    void run() override;
    void runOnIndices(const unsigned int* indices, unsigned int size_indices, CollisionTemplate collisionTemplate,
                      CudaStreamIndex streamIndex) override;

private:
    K17CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level);
};

#endif
