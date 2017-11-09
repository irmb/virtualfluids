/****************************************************************************
** $Id: qt/composer.cpp   3.3.1   edited May 27 2003 $
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#include "composer.h"
#include "smtp.h"

#include <qlineedit.h>
#include <qmultilineedit.h>
#include <qpushbutton.h>
#include <qlabel.h>
#include <qlayout.h>

Composer::Composer( QWidget *parent ) : QWidget( parent )
{
    outStandigEmails = 0;
    QGridLayout * layout = new QGridLayout( this, 1, 1, 6 );

    layout->addWidget( new QLabel( tr( "From:" ), this ), 0, 0 );
    from = new QLineEdit( this );
    layout->addWidget( from, 0, 1 );

    layout->addWidget( new QLabel( tr( "To:" ), this ), 1, 0 );
    to = new QLineEdit( this );
    layout->addWidget( to, 1, 1 );

    layout->addWidget( new QLabel( tr( "Subject:" ), this ), 2, 0 );
    subject = new QLineEdit( this );
    layout->addWidget( subject, 2, 1 );

    message = new QMultiLineEdit( this );
    layout->addMultiCellWidget( message, 3, 3, 0, 1 );

    send = new QPushButton( tr( "&Send" ), this );
    layout->addWidget( send, 4, 0 );
    connect( send, SIGNAL( clicked() ), this, SLOT( sendMessage() ) );

    sendStatus = new QLabel( this );
    layout->addWidget( sendStatus,4, 1 );
}
void Composer::sendMessage(const QString &from, const QString &to,const QString &subject,const QString &body)
{
   outStandigEmails++;

   Smtp *smtp = new Smtp( from,to, subject,body );
   connect( smtp, SIGNAL(status(const QString &)),this, SLOT(coutText(const QString &)) );
   connect( smtp, SIGNAL(destroyed()),this, SLOT(reduceEmailCounter()) );
}
void Composer::sendMessage(const QString &from, const QString &to,const QString &subject,const QString &body,const  QString &smtpServer)
{
   outStandigEmails++;

   Smtp *smtp = new Smtp( from,to, subject,body,smtpServer );
   connect( smtp, SIGNAL(status(const QString &)),this, SLOT(coutText(const QString &)) );
   connect( smtp, SIGNAL(destroyed()),this, SLOT(reduceEmailCounter()) );
}

void Composer::sendMessage()
{
   send->setEnabled( FALSE );
   sendStatus->setText( tr( "Looking up mail servers" ) );

   outStandigEmails++;
   Smtp *smtp = new Smtp(from->text(), to->text(),subject->text(),message->text());
   connect( smtp, SIGNAL(status(const QString &)),sendStatus, SLOT(setText(const QString &)) );
   connect( smtp, SIGNAL(destroyed()),this, SLOT(reduceEmailCounter()) );
   connect( smtp, SIGNAL(destroyed()),this, SLOT(enableSend()) );
}

void Composer::enableSend()
{
    send->setEnabled( TRUE );
}

void Composer::coutText(const QString& text)
{
   cout<<text<<endl;
}

void Composer::reduceEmailCounter()
{
    outStandigEmails--;
    cout<<"outstanding Emails : "<<outStandigEmails<<endl;
}
