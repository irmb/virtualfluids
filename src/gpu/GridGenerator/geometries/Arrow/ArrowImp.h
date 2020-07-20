#ifndef ArrowImp_H
#define ArrowImp_H

#include <memory>

#include "global.h"

#include "Arrow.h"

struct Vertex;

class ArrowImp : public Arrow 
{
public:
	VIRTUALFLUIDS_GPU_EXPORT virtual ~ArrowImp();
	VIRTUALFLUIDS_GPU_EXPORT static std::shared_ptr<Arrow> make(const Vertex &start, const Vertex &end);

	VIRTUALFLUIDS_GPU_EXPORT std::shared_ptr<Vertex> getStart() const;
	VIRTUALFLUIDS_GPU_EXPORT std::shared_ptr<Vertex> getEnd() const;

	VIRTUALFLUIDS_GPU_EXPORT void print() const;
private:
	ArrowImp(const Vertex &start, const Vertex &end);

	std::shared_ptr<Vertex> start;
	std::shared_ptr<Vertex> end;
};



#endif
