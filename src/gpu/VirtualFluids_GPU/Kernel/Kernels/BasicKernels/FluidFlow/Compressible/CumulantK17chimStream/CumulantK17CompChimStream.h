#ifndef CUMULANT_K17_COMP_CHIM_SPARSE_H
#define CUMULANT_K17_COMP_CHIM_SPARSE_H

#include "Kernel/KernelImp.h"

class CumulantK17CompChimStream : public KernelImp
{
public:
    static std::shared_ptr<CumulantK17CompChimStream> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run() override;
    void runOnIndices(const unsigned int *indices, unsigned int size_indices, CudaStreamIndex streamIndex=CudaStreamIndex::Legacy) override;

private:
    CumulantK17CompChimStream();
    CumulantK17CompChimStream(std::shared_ptr<Parameter> para, int level);
};

#endif 
