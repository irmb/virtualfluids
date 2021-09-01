#ifndef CUMULANT_K17_COMP_CHIM_SPARSE_H
#define CUMULANT_K17_COMP_CHIM_SPARSE_H

#include "Kernel/KernelImp.h"

class CumulantK17CompChimSparse : public KernelImp
{
public:
    static std::shared_ptr<CumulantK17CompChimSparse> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run() override;
    void runOnIndices(const unsigned int *indices, unsigned int size_indices, int stream = -1) override;

private:
    CumulantK17CompChimSparse();
    CumulantK17CompChimSparse(std::shared_ptr<Parameter> para, int level);
};

#endif 
