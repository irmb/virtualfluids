//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBTABLEMODEL_H
#define UBTABLEMODEL_H

#include <iostream>

class UbTableModel 
{
public:
   const static int COL_TYPE_STRING  = 1;
   const static int COL_TYPE_BOOL    = 2;
   const static int COL_TYPE_INTEGER = 3;
   const static int COL_TYPE_DOUBLE  = 4;

	UbTableModel();
	virtual ~UbTableModel();

	//////////////////////////////////////////////////////////////////////////
	//void objectChanged(UbObservable*);
	//void objectWillBeDeleted(UbObservable*);

	virtual int getColumnNumber() = 0;
	virtual int getRowNumber()    = 0;
   virtual std::string getColumnLabel(int column) = 0;
	virtual int getColumnType(int column) = 0;
	virtual std::string getStringValue(int row, int col) = 0;
   virtual void setStringValue(int row, int col, std::string str) = 0;
	virtual int getSelectedRowIndex() = 0;
	//virtual bool GetBoolValue(int row, int col) = 0;
};

#endif //UBTABLEMODEL_H
