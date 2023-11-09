#ifndef K15CompressibleNavierStokesPorousMedia_H
#define K15CompressibleNavierStokesPorousMedia_H

#include "Kernel/KernelImp.h"

class PorousMedia;

class K15CompressibleNavierStokesPorousMedia : public KernelImp
{
public:
	static std::shared_ptr<K15CompressibleNavierStokesPorousMedia> getNewInstance(std::shared_ptr< Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level);
	void run();

private:
	K15CompressibleNavierStokesPorousMedia();
	K15CompressibleNavierStokesPorousMedia(std::shared_ptr< Parameter> para, std::vector<std::shared_ptr<PorousMedia>> pm, int level);

	std::vector<std::shared_ptr<PorousMedia>> pm;
};

#endif 