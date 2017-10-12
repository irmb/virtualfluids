#ifndef ArrowImp_H
#define ArrowImp_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <memory>
#include "Arrow.h"

struct Vertex;

class ArrowImp : public Arrow 
{
public:
	GridGenerator_EXPORT virtual ~ArrowImp();
	GridGenerator_EXPORT static std::shared_ptr<Arrow> make(const Vertex &start, const Vertex &end);

	GridGenerator_EXPORT std::shared_ptr<Vertex> getStart() const;
	GridGenerator_EXPORT std::shared_ptr<Vertex> getEnd() const;

	GridGenerator_EXPORT void print() const;
private:
	ArrowImp(const Vertex &start, const Vertex &end);

	std::shared_ptr<Vertex> start;
	std::shared_ptr<Vertex> end;
};



#endif
