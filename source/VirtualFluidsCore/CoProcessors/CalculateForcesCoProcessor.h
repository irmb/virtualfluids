/*
 *  D3Q27ForcesCoProcessor.h
 *
 *  Created on: 29.09.2012
 *  Author: K. Kucher
 */

#ifndef D3Q27ForcesCoProcessor_H
#define D3Q27ForcesCoProcessor_H

#include <memory>
#include <string>
#include <vector>

#include "CoProcessor.h"

class ForceCalculator;
class Communicator;
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;

class CalculateForcesCoProcessor: public CoProcessor 
{
public:
   //! Constructor
   //! \param v - velocity of fluid in LB units
   //! \param a - area of object in LB units
   CalculateForcesCoProcessor(std::shared_ptr<Grid3D> grid, std::shared_ptr<UbScheduler> s,
                            const std::string &path,
       std::shared_ptr<Communicator> comm, double v, double a);
	virtual ~CalculateForcesCoProcessor();             
	void process(double step); 
   void addInteractor(std::shared_ptr<D3Q27Interactor> interactor);
protected:
	void collectData(double step);
   void calculateForces();
   UbTupleDouble3 getForces(int x1, int x2, int x3, std::shared_ptr<DistributionArray3D> distributions, std::shared_ptr<BoundaryConditions> bc);
   void calculateCoefficients();
   void write(std::ofstream *fileObject, double value, char *separator);
private:
   std::string path;
   std::shared_ptr<Communicator> comm;
   std::vector<std::shared_ptr<D3Q27Interactor> > interactors;
   double forceX1global;
   double forceX2global;
   double forceX3global;
   double v;     //!< is the speed of the object relative to the fluid
   double a;     //!< is the reference area
   double C1;
   double C2;
   double C3;

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive & ar, const unsigned int version)
   {
      ar & boost::serialization::base_object<CoProcessor>(*this);
      ar & path;
      ar & v;
      ar & a;
      ar & interactors;
   }
};


#endif /* D3Q27ForcesCoProcessor_H */
