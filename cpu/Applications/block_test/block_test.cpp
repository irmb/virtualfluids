#include "block_test_incompressible.hpp"

int main(int argc, char* argv[])
{

   if ( argv != NULL )
   {
      if (argc > 1)
      {
         block_test_incompressible(argv[1], argv[2]);
      }
      else
      {
         cout << "Configuration file must be set!: " <<  argv[0] << " <config file>" << endl << std::flush;
      }
   }

   return 0;
}

