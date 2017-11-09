#include <qapplication.h>
#include <qsocket.h>
#include <qstring.h>

#include "composer.h"
#include "smtp.h"
#include <iostream>
using namespace std;


int main( int argc, char **argv )
{
    QApplication a( argc, argv );
    
    
    Composer *c = new Composer();
    //a.setMainWidget(c);
    //c->show();
    //return a.exec();

    QString from    = "muffmolch@gmx.de";
    QString to      = "muffmolch@gmx.de";
    QString subject = "Kaffeeautomat-Reinigung die ";
    QString body    = "Danke für den Gratis-Kaffee";
    
    int j=1;
    for(int i=0;i<2;i++)
    {
       QString subject2;
       subject2= subject+QString::number(j++)+".";
       c->sendMessage(from,to,subject2,body,"mail.irmb.bau.tu-bs.de");
    }

    return a.exec();
}
