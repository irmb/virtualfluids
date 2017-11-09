//#define MPI_LOGGING


#include <mpi.h>
#if defined( MPI_LOGGING )
	#include <mpe.h>
#endif

#include <string>
#include <iostream>

#include "VirtualFluids_GPU/Application/App.h"
using namespace std;

int main( int argc, char* argv[]){
   MPI_Init(&argc, &argv);
   string str, str2; 
   if ( argv != NULL )
   {
      str = static_cast<string>(argv[0]);
      if (argc > 1)
      {
         str2 = static_cast<string>(argv[1]);
         App::getInstanz()->run(str2);
      }
      else
      {
         cout << "Configuration file must be set!: lbmgm <config file>" << endl << std::flush;
         //MPI_Abort(MPI_COMM_WORLD, -1);
      }
   }
   /*
   MPE_Init_log() & MPE_Finish_log() are NOT needed when
   liblmpe.a is linked with this program.  In that case,
   MPI_Init() would have called MPE_Init_log() already.
   */
#if defined( MPI_LOGGING )
   MPE_Init_log();
#endif

   

#if defined( MPI_LOGGING )
   if ( argv != NULL )
      MPE_Finish_log( argv[0] );
   if ( str != "" )
      MPE_Finish_log( str.c_str() );
   else
      MPE_Finish_log( "TestLog" );
#endif

   MPI_Finalize();
   return 0;
}
