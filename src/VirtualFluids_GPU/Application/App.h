#ifndef _APP_H_
#define _APP_H_
#include <string>

#include <VirtualFluidsDefinitions.h>

////Gemeinsame Variablen
//Buffer2D <float> sbuf_t; 
//Buffer2D <float> rbuf_t;
//Buffer2D <float> sbuf_b;
//Buffer2D <float> rbuf_b;
//
////for MPI
//Communicator* comm;
//int numprocs, myid;
//
//WriteLog output;

class App
{
public:
	VF_PUBLIC static App* getInstanz();
    VF_PUBLIC void run(std::string &cstr);
protected:

private:
   static App* instanz;
   App(){}
   App(const App&);
};

#endif
