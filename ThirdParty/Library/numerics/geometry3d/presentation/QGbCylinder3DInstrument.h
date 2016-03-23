#ifndef QGBCYLINDER3DINSTRUMENT_H
#define QGBCYLINDER3DINSTRUMENT_H

#include <QtGui/QDialog>
#include <QtGui/QWidget>

#include "./QGbCylinder3DInstrumentUI.h"
#include "./QGbObject3DInstrument.h"

class GbCylinder3D;
class GbObject3D;

class QGbCylinder3DInstrument : public QDialog
{

   Q_OBJECT

public:
   QGbCylinder3DInstrument( QWidget* parent = 0, Qt::WFlags flags = 0 );
   ~QGbCylinder3DInstrument();
   void setGbCylinder3D(GbCylinder3D* cylinder);
   GbCylinder3D* getGbCylinder3D();

protected:
   GbCylinder3D* gbCylinder;

private:
   Ui::QGbCylinder3DInstrument ui;

private slots:
   void on_pBtnOK_clicked();
   void on_pBtnCancel_clicked();
};

#endif