#include "smtp.h"

#include <qtextstream.h>
#include <qsocket.h>
#include <qdns.h>
#include <qtimer.h>
#include <qapplication.h>
#include <qmessagebox.h>
#include <qregexp.h>

#include <iostream>
using namespace std;

Smtp::Smtp( const QString &from, const QString &to,const QString &subject,const QString &body )
{
   mSocket = new QSocket(this); 
   connect ( mSocket, SIGNAL( readyRead() ),this, SLOT( readyRead() ) );
   connect ( mSocket, SIGNAL( connected() ),this, SLOT( connected() ) );
   
   mTextStream = new QTextStream(mSocket);

   //wenn man es nicht über autom. ermittlung machen möchte
   mMxLookup = new QDns( from.mid( from.find( '@' )+1 ), QDns::Mx );
   connect( mMxLookup, SIGNAL(resultsReady()),this, SLOT(dnsLookupHelper()) );

   mFrom    = from;
   mTo      = to;
   mSubject = subject;
   mBody    = body;

   state = Init;
}

Smtp::Smtp( const QString &from, const QString &to,const QString &subject,const QString &body,const  QString &smtpServer )
{
   mSocket = new QSocket(this); 
   connect ( mSocket, SIGNAL( readyRead() ),this, SLOT( readyRead() ) );
   connect ( mSocket, SIGNAL( connected() ),this, SLOT( connected() ) );
   
   mTextStream = new QTextStream(mSocket);
   emit status( tr( "Connecting to " + smtpServer ) );
   mSocket->connectToHost(smtpServer,25);

   mFrom    = from;
   mTo      = to;
   mSubject = subject;
   mBody    = body;

   state = Init;
}

Smtp::~Smtp()
{
   delete mTextStream;
   delete mSocket;
}

void Smtp::dnsLookupHelper()
{
   QValueList<QDns::MailServer> s = mMxLookup->mailServers();
   
   if ( s.isEmpty() ) 
   {
      if ( !mMxLookup->isWorking() )  emit status( tr( "Smtp::dnsLookupHelper() Error in MX record lookup" ) );
      return;
   }

   emit status( tr( "Connecting to %1" ).arg( s.first().name ) );
   mSocket->connectToHost( s.first().name, 25 );
}


void Smtp::connected()
{
   emit status( tr( "Connected to %1" ).arg( mSocket->peerName() ) );
}

void Smtp::readyRead()
{
   // SMTP is line-oriented
   if ( !mSocket->canReadLine() )  return;

   QString responseLine;
   QString response;

   do 
   {
      responseLine = mSocket->readLine();
      response += responseLine;
   } while( mSocket->canReadLine() && responseLine[3] != ' ' );
   
   responseLine.truncate( 3 );
   
   if ( state == Init && responseLine[0] == '2' ) 
   {
      // banner was okay, let's go on
      *mTextStream << "HELO there\r\n";
      state = Mail;
   } 
   else if ( state == Mail && responseLine[0] == '2' ) 
   {
      // HELO response was okay (well, it has to be)
      *mTextStream<<"MAIL FROM: <"<<mFrom<<">\r\n";
      state = Rcpt;
   } 
   else if ( state == Rcpt && responseLine[0] == '2' ) 
   {
      *mTextStream << "RCPT TO: <" << mTo << ">\r\n";
      state = Data;
   }
   else if ( state == Data && responseLine[0] == '2' ) 
   {
      *mTextStream << "DATA\r\n";
      state = Body;
   } 
   else if ( state == Body && responseLine[0] == '3' )
   {
      QString message = QString::fromLatin1( "From: " ) + mFrom +
                        QString::fromLatin1( "\nTo: " ) + mTo +
                        QString::fromLatin1( "\nSubject: " ) + mSubject +
                        QString::fromLatin1( "\n\n" ) + mBody + "\n";
      message.replace( QString::fromLatin1( "\n" ),QString::fromLatin1( "\r\n" ) );
      message.replace( QString::fromLatin1( "\r\n.\r\n" ),QString::fromLatin1( "\r\n..\r\n" ) );
      
      *mTextStream << message << ".\r\n";
      state = Quit;
   } 
   else if ( state == Quit && responseLine[0] == '2' ) 
   {
      *mTextStream << "QUIT\r\n";
      // here, we just close.
      state = Close;
      emit status( tr( "Smtp::readyRead()::Message sent" ) );
   } 
   else if ( state == Close ) 
   {
      deleteLater();
      return;
   } 
   else 
   {
      // something broke.
      QMessageBox::warning( qApp->activeWindow(),tr( "Qt Mail Example" ),
         tr( "Unexpected reply from SMTP server:\n\n" ) +response );

      state = Close;
   }
}
