#ifndef CUMULANT_K17_ALMIGHTY_H
#define CUMULANT_K17_ALMIGHTY_H

#include "Kernel/KernelImp.h"
#include "Parameter/Parameter.h"

template<TurbulenceModel turbulenceModel> 
class CumulantK17Almighty : public KernelImp
{
public:
	static std::shared_ptr< CumulantK17Almighty<turbulenceModel> > getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run() override;
    void runOnIndices(const unsigned int *indices, unsigned int size_indices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex) override;

private:
    CumulantK17Almighty();
    CumulantK17Almighty(std::shared_ptr<Parameter> para, int level);
};

#endif 
