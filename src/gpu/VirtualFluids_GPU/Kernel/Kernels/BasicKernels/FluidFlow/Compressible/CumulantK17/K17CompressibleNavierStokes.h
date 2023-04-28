#ifndef CUMULANT_K17_H
#define CUMULANT_K17_H

#include "Kernel/KernelImp.h"
#include "Parameter/Parameter.h"

template<TurbulenceModel turbulenceModel> 
class K17CompressibleNavierStokes : public KernelImp
{
public:
	static std::shared_ptr< K17CompressibleNavierStokes<turbulenceModel> > getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run() override;
    void runOnIndices(const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex) override;

private:
    K17CompressibleNavierStokes();
    K17CompressibleNavierStokes(std::shared_ptr<Parameter> para, int level);
};

#endif 
