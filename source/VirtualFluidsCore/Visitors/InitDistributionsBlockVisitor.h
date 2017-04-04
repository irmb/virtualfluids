#ifndef InitDistributionsBlockVisitor_H
#define InitDistributionsBlockVisitor_H

#include <boost/foreach.hpp>

#include "Block3DVisitor.h"
#include "D3Q27System.h"
#include "Block3D.h"

#include <MuParser/include/muParser.h>

/*================================================================================*/
/*  D3Q27ETInitDistributionsBlockVisitor                                             */
/*                                                                                */
/**
more flexible way to initialize flow area
you can define functions to calculate macroscopic values for feq 
!!! x1,x2,x3 are automatically defined via this adapter and are the real world
vertex coordinates !!!

if function is invalid an UbException with detailed information is thrown

<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 19.04.08
*/ 

//! \details example:<BR>
//! D3Q27InitDistributionsBlockVisitor init(1.0,0.0,0.0,0.0);<BR>
//! Bem.: rho=0.0 bei inkompressibel<BR>
//! init.setVx1("0.01*x2");<BR>
//! init.setVx1("0.01*x2^2");<BR>
//! patch.adaptByPatchCriterion(init);



class InitDistributionsBlockVisitor : public Block3DVisitor
{
public:
   typedef std::numeric_limits<LBMReal> D3Q27RealLim;

public:
   InitDistributionsBlockVisitor();
   //D3Q27ETInitDistributionsBlockVisitor(LBMReal rho, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
   //! Constructor
   //! \param nu - viscosity
   //! \param rho - density
   //! \param vx1 - velocity in x
   //! \param vx2 - velocity in y
   //! \param vx3 - velocity in z
   InitDistributionsBlockVisitor( LBMReal nu, LBMReal rho, LBMReal vx1=0.0, LBMReal vx2=0.0, LBMReal vx3=0.0);
   //////////////////////////////////////////////////////////////////////////
   //automatic vars are: x1,x2, x3
   //ussage example: setVx1("x1*0.01+x2*0.003")
   //////////////////////////////////////////////////////////////////////////
   void setVx1( const mu::Parser& parser);
   void setVx2( const mu::Parser& parser);
   void setVx3( const mu::Parser& parser);
   void setRho( const mu::Parser& parser);

   void setVx1( const std::string& muParserString);
   void setVx2( const std::string& muParserString);
   void setVx3( const std::string& muParserString);
   void setRho( const std::string& muParserString);
   //////////////////////////////////////////////////////////////////////////
   void setVx1( LBMReal vx1 );
   void setVx2( LBMReal vx2 );
   void setVx3( LBMReal vx3 );
   void setRho( LBMReal rho );
   void setNu( LBMReal nu );

   void visit(Grid3DPtr grid, Block3DPtr block);

protected:
   void checkFunction(mu::Parser fct);

private:
   mu::Parser muVx1;
   mu::Parser muVx2;
   mu::Parser muVx3;
   mu::Parser muRho;
   LBMReal nu;
};

#endif //D3Q27INITDISTRIBUTIONSPATCHVISITOR_H
