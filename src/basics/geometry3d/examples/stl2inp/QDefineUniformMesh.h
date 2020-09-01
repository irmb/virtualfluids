#ifndef QDEFINEUNIFORMMESH_H
#define QDEFINEUNIFORMMESH_H

#include <QtGui/QDialog>
#include "./QDefineUniformMeshUI.h"

class QDefineUniformMesh : public QDialog
{
	Q_OBJECT

public:
	QDefineUniformMesh(QWidget *parent = 0, Qt::WFlags flags = 0);
	~QDefineUniformMesh();

	void	setStartLevel(int startLevel)	{ ui.spinBox_startLevel->setValue(startLevel); }
	void	setStopLevel(int stopLevel)		{ ui.spinBox_stopLevel->setValue(stopLevel); }
	void	setDelta(double delta)		{ ui.doubleSpinBox_delta->setValue(delta); }
	void	setNX1(int nx1)				{ ui.spinBox_nx1->setValue(nx1); }
	void	setNX2(int nx2)				{ ui.spinBox_nx2->setValue(nx2); }
	void	setNX3(int nx3)				{ ui.spinBox_nx3->setValue(nx3); }

	int		getStartLevel()	{ return ui.spinBox_startLevel->value(); }
	int		getStopLevel()	{ return ui.spinBox_stopLevel->value(); }
	double	getDelta()		{ return ui.doubleSpinBox_delta->value(); }
	int		getNX1()		{ return ui.spinBox_nx1->value(); }
	int		getNX2()		{ return ui.spinBox_nx2->value(); }
	int		getNX3()		{ return ui.spinBox_nx3->value(); }

private:
	Ui::QDefineUniformMesh ui;

//private slots:
};

#endif // QDEFINEUNIFORMMESH_H
