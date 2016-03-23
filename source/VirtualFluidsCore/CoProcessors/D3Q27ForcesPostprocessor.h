/*
 *  D3Q27ForcesPostprocessor.h
 *
 *  Created on: 29.09.2012
 *  Author: K. Kucher
 */

#ifndef D3Q27ForcesPostprocessor_H
#define D3Q27ForcesPostprocessor_H

#include "Postprocessor.h"
#include "Communicator.h"
#include "D3Q27Interactor.h"

class D3Q27ForcesPostprocessor: public Postprocessor 
{
public:
   //! Constructor
   //! \param v - velocity of fluid in LB units
   //! \param a - area of object in LB units
   D3Q27ForcesPostprocessor(Grid3DPtr grid, UbSchedulerPtr s,
                            const std::string &path,
                            CommunicatorPtr comm, double v, double a);
	virtual ~D3Q27ForcesPostprocessor();             
	void update(double step); 
   void addInteractor(D3Q27InteractorPtr interactor);
protected:
	void collectPostprocessData(double step);
   void calculateForces();
   UbTupleDouble3 getForces(int x1, int x2, int x3, DistributionArray3DPtr distributions, D3Q27BoundaryConditionPtr bc);
   void calculateCoefficients();
   void write(std::ofstream *fileObject, double value, char *separator);
private:
   std::string path;
   CommunicatorPtr comm;
   std::vector<D3Q27InteractorPtr> interactors;
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
      ar & boost::serialization::base_object<Postprocessor>(*this);
      ar & path;
      ar & v;
      ar & a;
      ar & interactors;
   }
};


#endif /* D3Q27ForcesPostprocessor_H */
