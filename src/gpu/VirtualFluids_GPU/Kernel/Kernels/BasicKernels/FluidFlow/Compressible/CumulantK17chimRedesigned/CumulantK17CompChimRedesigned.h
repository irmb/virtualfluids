#ifndef CUMULANT_K17_COMP_CHIM_REDESIGN_H
#define CUMULANT_K17_COMP_CHIM_REDESIGN_H

#include "Kernel/KernelImp.h"

class CumulantK17CompChimRedesigned : public KernelImp
{
public:
    static std::shared_ptr<CumulantK17CompChimRedesigned> getNewInstance(std::shared_ptr<Parameter> para, int level);
	void run() override;
    void runOnIndices(const unsigned int *indices, unsigned int size_indices, int stream = -1) override;

private:
    CumulantK17CompChimRedesigned();
    CumulantK17CompChimRedesigned(std::shared_ptr<Parameter> para, int level);
};

#endif 
