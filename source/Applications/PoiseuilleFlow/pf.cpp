#include <iostream>
#include <string>

#include "pf.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
   try
   {
      pf1();
      return 0;
   }
   catch (std::exception& e)
   {
      cerr << e.what() << endl << flush;
   }
   catch (std::string& s)
   {
      cerr << s << endl;
   }
   catch (...)
   {
      cerr << "unknown exception" << endl;
   }
}
