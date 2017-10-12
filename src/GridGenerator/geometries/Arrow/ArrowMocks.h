#ifndef ArrowMocks_H
#define ArrowMocks_H

#include "GridGenerator/global.h"
#include <memory>

#include "Arrow.h"
#include "../Vertex/Vertex.cuh"

class ArrowStub : public Arrow
{
public:
	virtual ~ArrowStub() {}
	static std::shared_ptr<Arrow> make(std::shared_ptr<Vertex> start, std::shared_ptr<Vertex> end) 
	{
		return std::shared_ptr<Arrow>(new ArrowStub(start, end));
	}

	std::shared_ptr<Vertex> getStart() const { return this->start; }
	std::shared_ptr<Vertex> getEnd() const { return this->end; };


	virtual void print() const override
	{
		throw std::logic_error("The method or operation is not implemented.");
	}

private:
	ArrowStub(std::shared_ptr<Vertex> start, std::shared_ptr<Vertex> end) : start(start), end(end) {}

	std::shared_ptr<Vertex> start;
	std::shared_ptr<Vertex> end;
};






#endif
