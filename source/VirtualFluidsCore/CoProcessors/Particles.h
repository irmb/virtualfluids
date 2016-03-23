#ifndef PARTICLES_H
#define PARTICLES_H

#include <boost/shared_ptr.hpp>
# define numOfvarClass 15

class Particles;
typedef boost::shared_ptr<Particles> ParticlesPtr;

class Particles 
{
public:
	double x,y,z,xold,yold,zold;
	double vxold,vyold,vzold;
	double vxoldf,vyoldf,vzoldf;
	int ID;
	int rankOfParticle;
    int level;
	//double dx;


	Particles();
	~Particles();
protected:
private:
	
};
#endif

