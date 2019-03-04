#ifndef BGK_PLUS_COMP_SP27_H
#define BGK_PLUS_COMP_SP27_H

#include "Kernel\Kernel.h"

#include <memory>

class BGKPlusCompSP27 : public Kernel
{
public:
	static std::shared_ptr<Kernel> getNewInstance(std::shared_ptr< Parameter> para, int level);
	void run();
	bool checkParameter();



private:
	BGKPlusCompSP27();
	BGKPlusCompSP27(std::shared_ptr< Parameter> para, int level);
	std::shared_ptr< Parameter> para;
	int level;
};

#endif 