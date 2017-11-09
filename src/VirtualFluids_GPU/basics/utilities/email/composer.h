#ifndef COMPOSER_H
#define COMPOSER_H

#include <iostream>
using namespace std;

#include <qwidget.h>

class QLineEdit;
class QMultiLineEdit;
class QLabel;
class QPushButton;


class Composer : public QWidget
{
    Q_OBJECT
public:
    Composer( QWidget *parent = 0 );
    void sendMessage(const QString &from, const QString &to,const QString &subject,const QString &body);
    void sendMessage(const QString &from, const QString &to,const QString &subject,const QString &body,const  QString &smtpServer);

private slots:
    void sendMessage();
    void enableSend();
    void coutText(const QString& text);
    void reduceEmailCounter();

private:
   int outStandigEmails; 
   QLineEdit *from, *to, *subject;
   QMultiLineEdit *message;
   QLabel * sendStatus;
   QPushButton * send;
};


#endif
