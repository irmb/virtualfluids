#include "DensityBCAdapter.h"
#include "basics/utilities/UbLogger.h"
#include "basics/utilities/UbInfinity.h"

using namespace std;
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const double& dens, const double& startTime, const double& endTime )
{
   this->densBCs.push_back( BCFunction(dens,startTime,endTime) );
   this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const BCFunction& densBC )
{
   this->densBCs.push_back(densBC);
   this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const std::vector< BCFunction >& densBCs)
{
   this->densBCs = densBCs;
   this->init();
}
/*==========================================================*/
DensityBCAdapter::DensityBCAdapter(const mu::Parser& function, const double& startTime, const double& endTime )
{
   this->densBCs.push_back(BCFunction(function,startTime,endTime));
   this->init();
}
/*==========================================================*/
void DensityBCAdapter::init()
{
   this->timeStep = 0.0;

   this->x1 = 0.0;
   this->x2 = 0.0;
   this->x3 = 0.0;

   this->tmpDensityFunction = NULL;

   try //initilialization and validation of functions
   {
      for(size_t pos=0; pos<densBCs.size(); ++pos)
      {
         if( !(    UbMath::equal( BCFunction::INFCONST, densBCs[pos].getEndTime() )
                && UbMath::greaterEqual( this->timeStep,  densBCs[pos].getStartTime()  ) ) )
         { 
            this->setTimeDependent();
         }

         densBCs[pos].getFunction().DefineVar("t" , &this->timeStep);
         densBCs[pos].getFunction().DefineVar("x1", &this->x1      );
         densBCs[pos].getFunction().DefineVar("x2", &this->x2      );
         densBCs[pos].getFunction().DefineVar("x3", &this->x3      );

         densBCs[pos].getFunction().Eval(); //<-- validation
      }
   }
   catch(mu::Parser::exception_type& e){ stringstream error; error<<"mu::parser exception occurs, message("<<e.GetMsg()<<"), formula("<<e.GetExpr()+"), token("+e.GetToken()<<")"
                                         <<", pos("<<e.GetPos()<<"), error code("<<e.GetCode(); throw UbException(error.str()); }
   catch(...)                          { throw UbException(UB_EXARGS,"unknown exception" );                       }
}
/*==========================================================*/
void DensityBCAdapter::init(const D3Q27Interactor* const& interactor, const double& time)
{
   this->timeStep           = time;
   this->tmpDensityFunction = NULL;
   double maxEndtime        = -Ub::inf;

   //aktuelle Densityfunction bestimmen
   for(size_t pos=0; pos<densBCs.size(); ++pos)
   {
      if( UbMath::equal(densBCs[pos].getEndTime(),BCFunction::INFTIMEDEPENDENT)) maxEndtime=Ub::inf;
      maxEndtime = UbMath::max(maxEndtime,densBCs[pos].getStartTime(),densBCs[pos].getEndTime()); //startTime abfragen, da  INFCONST=-10

      if( UbMath::greaterEqual(this->timeStep,densBCs[pos].getStartTime()) ) 
      {
         if(   UbMath::lessEqual(this->timeStep,densBCs[pos].getEndTime())
            || UbMath::equal(densBCs[pos].getEndTime(),(double)BCFunction::INFCONST)
            || UbMath::equal(densBCs[pos].getEndTime(),(double)BCFunction::INFTIMEDEPENDENT) )
         {
            tmpDensityFunction = &densBCs[pos].getFunction();
            break;
         }
      }
   }

   //wenn funktionen zweitlich konstant sind und bis t=unendlich gelten
   //kann man zeitabhaengigkeit deaktivieren
   if( UbMath::greaterEqual(time,maxEndtime) ) this->unsetTimeDependent();

   UBLOG(logDEBUG4,"D3Q27DensityBCAdapter::init(time="<<time<<") "
                    <<", rho= \""<<(tmpDensityFunction ? tmpDensityFunction->GetExpr() : "-")
                    <<"\", timedependant="<<(this->isTimeDependent() ? "true" : "false") );
}
/*==========================================================*/
void DensityBCAdapter::update( const D3Q27Interactor* const& interactor, const double& time ) 
{
   this->init(interactor,time);
}
/*==========================================================*/
void DensityBCAdapter::adaptBCForDirection( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& q, const int& fdirection, const double& time )
{
   bc->setDensityBoundaryFlag(D3Q27System::INVDIR[fdirection],secondaryBcOption);
   bc->setQ((float)q,fdirection);
}
/*==========================================================*/
void DensityBCAdapter::adaptBC( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& time )
{
   this->setNodeDensity(interactor,bc,worldX1,worldX2,worldX3,time);
   bc->setBcAlgorithmType(algorithmType);
}
/*==========================================================*/
void DensityBCAdapter::setNodeDensity( const D3Q27Interactor& interactor, SPtr<BoundaryConditions> bc, const double& worldX1, const double& worldX2, const double& worldX3, const double& timestep) 
{
   //Geschwindigkeiten setzen
   try
   {
      //PunktKoordinaten bestimmen
      this->x1 = worldX1;
      this->x2 = worldX2;
      this->x3 = worldX3;
      this->timeStep = timestep;

      if(tmpDensityFunction) bc->setBoundaryDensity((float)tmpDensityFunction->Eval());  
   }
   catch(mu::Parser::exception_type& e){ stringstream error; error<<"mu::parser exception occurs, message("<<e.GetMsg()<<"), formula("<<e.GetExpr()+"), token("+e.GetToken()<<")"
                                          <<", pos("<<e.GetPos()<<"), error code("<<e.GetCode(); throw UbException(error.str()); }
   catch(...)                          { throw UbException(UB_EXARGS,"unknown exception" ); }
}
/*==========================================================*/
double DensityBCAdapter::getDensity(const double& x1, const double& x2, const double& x3, const double& timeStep)
{
   this->x1 = x1;
   this->x2 = x2;
   this->x3 = x3;
   this->timeStep = timeStep;

   if(!tmpDensityFunction) return 0.0;

   return tmpDensityFunction->Eval();  
}
/*==========================================================*/
string DensityBCAdapter::toString()
{
   stringstream info;
   info<<"D3Q27DensityBCAdapter:\n";
   info<<" #dens-functions = "<<(int)densBCs.size()<<endl;
   info<<" protected variables: x1, x2, x3, t"<<endl;

   for(size_t i=0; i<densBCs.size(); ++i)
   {
      info<<"\n   dens-function nr."<<i<<":"<<endl;
      info<<densBCs[i].toString()<<endl;
   }
   return info.str();
}
