#ifndef ArrowImp_H
#define ArrowImp_H

#include <memory>

#include "global.h"
#include "GridGenerator_export.h"

#include "Arrow.h"

struct Vertex;

class ArrowImp : public Arrow 
{
public:
	GRIDGENERATOR_EXPORT virtual ~ArrowImp();
	GRIDGENERATOR_EXPORT static std::shared_ptr<Arrow> make(const Vertex &start, const Vertex &end);

	GRIDGENERATOR_EXPORT std::shared_ptr<Vertex> getStart() const;
	GRIDGENERATOR_EXPORT std::shared_ptr<Vertex> getEnd() const;

	GRIDGENERATOR_EXPORT void print() const;
private:
	ArrowImp(const Vertex &start, const Vertex &end);

	std::shared_ptr<Vertex> start;
	std::shared_ptr<Vertex> end;
};



#endif
