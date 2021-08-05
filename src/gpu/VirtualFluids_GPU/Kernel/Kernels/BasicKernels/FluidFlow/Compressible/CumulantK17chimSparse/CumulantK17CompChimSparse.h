#ifndef CUMULANT_K17_COMP_CHIM_SPARSE_H
#define CUMULANT_K17_COMP_CHIM_SPARSE_H

#include "Kernel/KernelImp.h"

class CumulantK17CompChimSparse : public KernelImp
{
public:
    static std::shared_ptr<CumulantK17CompChimSparse> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run();
    void runOnIndices(const unsigned int *indices, unsigned int size_indices) override;

private:
    CumulantK17CompChimSparse();
    CumulantK17CompChimSparse(std::shared_ptr<Parameter> para, int level);
    std::unique_ptr<std::pair<dim3, dim3>> calcGridDimensions(unsigned int size_Mat);
};

#endif 
