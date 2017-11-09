/****************************************************************************
** $Id: qt/smtp.h   3.3.1   edited May 27 2003 $
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#ifndef SMTP_H
#define SMTP_H

#include <qobject.h>
#include <qstring.h>

class QSocket;
class QTextStream;
class QDns;

class Smtp : public QObject
{
    Q_OBJECT

public:
    Smtp( const QString &from, const QString &to,const QString &subject, const QString &body );
    Smtp( const QString &from, const QString &to,const QString &subject, const QString &body,const  QString &smtpServer );
   ~Smtp();

signals:
    void status( const QString & );
    
private slots:
    void dnsLookupHelper();
    void readyRead();
    void connected();

private:
   enum State {Init,Mail,Rcpt,Data,Body,Quit,Close};
   int state;

   QString mFrom;
   QString mTo;
   QString mSubject;
   QString mBody;
   
   QSocket*     mSocket;
   QTextStream* mTextStream;
   QDns*        mMxLookup;
};

#endif
