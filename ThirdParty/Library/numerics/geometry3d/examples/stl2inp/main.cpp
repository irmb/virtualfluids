#include <iostream>
#include <cstdlib>

#include <QtGui/QApplication>
#include "./stl2inp.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    STL2INP w;

	w.show();
    a.connect(&a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()));
    return a.exec();
}
