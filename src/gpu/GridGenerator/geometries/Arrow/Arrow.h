#ifndef Arrow_H
#define Arrow_H

#include <memory>

struct Vertex;

class Arrow 
{
public:
	virtual ~Arrow() {};
protected:
	Arrow() {};

public:
	virtual std::shared_ptr<Vertex> getStart() const = 0;
	virtual std::shared_ptr<Vertex> getEnd() const = 0;
	virtual void print() const = 0;
};



#endif
