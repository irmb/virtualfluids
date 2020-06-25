#ifndef ArrowImp_H
#define ArrowImp_H

#include <memory>

#include "GridGenerator/global.h"

#include "Arrow.h"

struct Vertex;

class ArrowImp : public Arrow 
{
public:
	VF_PUBLIC virtual ~ArrowImp();
	VF_PUBLIC static std::shared_ptr<Arrow> make(const Vertex &start, const Vertex &end);

	VF_PUBLIC std::shared_ptr<Vertex> getStart() const;
	VF_PUBLIC std::shared_ptr<Vertex> getEnd() const;

	VF_PUBLIC void print() const;
private:
	ArrowImp(const Vertex &start, const Vertex &end);

	std::shared_ptr<Vertex> start;
	std::shared_ptr<Vertex> end;
};



#endif
