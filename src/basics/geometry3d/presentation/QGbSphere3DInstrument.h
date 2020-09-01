#ifndef QGBSPHERE3DINSTRUMENT_H
#define QGBSPHERE3DINSTRUMENT_H

#include <QtGui/QDialog>
#include <QtGui/QWidget>

#include "./QGbSphere3DInstrumentUI.h"
#include "./QGbObject3DInstrument.h"

class GbSphere3D;
class GbObject3D;


class QGbSphere3DInstrument : public QDialog
{
	Q_OBJECT

public:
	QGbSphere3DInstrument( QWidget* parent = 0, Qt::WFlags flags = 0 );
	~QGbSphere3DInstrument();
	void setGbSphere3D(GbSphere3D* sphere);     
	GbSphere3D* getGbSphere3D();
	//void SetGbObject3D(GbObject3D*);

protected:
	GbSphere3D* gbSphere;

private:
	Ui::QGbSphere3DInstrument ui;

private slots:
	void on_pBtnOK_clicked();
	void on_pBtnCancel_clicked();
};

#endif   
