#ifndef QGBOBJECT3DINSTRUMENT_H
#define QGBOBJECT3DINSTRUMENT_H

#include <QDialog>
#include "./QGbObject3DInstrumentUI.h"


class GbObject3D;

class QGbObject3DInstrument : public QDialog
{
	Q_OBJECT

public:
	QGbObject3DInstrument( QWidget* parent = 0, Qt::WFlags flags = 0 );
	~QGbObject3DInstrument();
	void setGbObject3D(GbObject3D* gbObject);           
	GbObject3D* getGbObject3D();

protected:
	GbObject3D *gbObject3D;

private:
	Ui::QGbObject3DInstrument ui;

private slots:
	void on_pBtnOK_clicked();
	void on_pBtnCancel_clicked();
};
#endif   
