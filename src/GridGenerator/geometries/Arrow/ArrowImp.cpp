#include "ArrowImp.h"

#include "../Vertex/Vertex.h"



 std::shared_ptr<Arrow> ArrowImp::make(const Vertex &start, const Vertex &end)
{
	return std::shared_ptr<ArrowImp>(new ArrowImp(start, end));
}

ArrowImp::ArrowImp(const Vertex &start, const Vertex &end) : start(std::make_shared<Vertex>(start)), end(std::make_shared<Vertex>(end))
{

}

ArrowImp::~ArrowImp()
{

}

std::shared_ptr<Vertex> ArrowImp::getStart() const
{
	return this->start;
}

std::shared_ptr<Vertex> ArrowImp::getEnd() const
{
	return this->end;
}

void ArrowImp::print() const
{
	printf("v1: ");
	start->print();
	printf("v2: ");
	end->print();
}

