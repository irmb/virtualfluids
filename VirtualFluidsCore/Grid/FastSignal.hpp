#include <list>
#include <boost/function.hpp>

namespace fast_signal
{
   template <class Signature>
   struct signal: public std::list<boost::function<Signature> >{};
}

// Abkürzung für die Schleife zum Aufrufen
#define CALL_SIGNALS( SIG, REF, ARGS )\
   for( fast_signal::signal<SIG>::iterator it = REF.begin();\
   it != REF.end(); ++it )\
   (*it)ARGS;


//Eigentlich war das schon alles was wir brauchen um unsere sehr einfache Signal und Slot Implementation zu verwenden. 
//Die Syntax ist zwar nicht ganz gleich wie bei Boost.Signals aber das sollten wir für die Zeitersparnis verkraften können. 
//Jetzt also ein einfaches Verwendungsbeispiel:
//#include <iostream>
//// ... Code von vorher
//
//// Testfunktion zum Registrieren als Slots
//void test( int& i )
//{
//   i += 1;
//}
//
//int main()
//{
//   tp_fast_signal::signal<void(int&)> test_signal; // Signal erstellen
//
//   // 2 Slots registrieren:
//   test_signal.push_back(&test);
//   test_signal.push_back(&test);
//
//   int test_int = 0;
//
//   // Alle Slots aufrufen
//   CALL_SIGNALS( void(int&), test_signal, (test_int) );
//
//   // Und noch testweise den aktuellen Wert von test_int ausgeben
//   std::cout << "test_int=" << test_int << std::endl;
//
//   return 0;
//}
