#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>
#include <math.h>

// Get current date/time, format is DD.MM.YYYY/HH:mm:ss
const std::string currentDateTime() {
   time_t     now = time(0);
   struct tm  tstruct;
   char       buf[80];
   tstruct = *localtime(&now);
   // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
   // for more information about date/time format
   strftime(buf, sizeof(buf), "%d.%m.%Y/%X", &tstruct);

   return buf;
}
//////////////////////////////////////////////////////////////////////////
void print(const std::string& s)
{
   std::cout << currentDateTime() << ": " << s <<std::endl<<std::flush;
}
//////////////////////////////////////////////////////////////////////////
//convert from double to int
//alles falsch hier
int cint(double x)
{
   double intpart;
   if (modf(x,&intpart)>=.5)
      return (int)floor(x)+1;
   else
   {
      if (x < 0)
         return (int)ceil(x);
      else
         return (int)floor(x);
   }
}

