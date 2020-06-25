#include <basics/utilities/UbTableModel.h>

UbTableModel::UbTableModel()
{
}

UbTableModel::~UbTableModel()
{
	//this->notifyObserversObjectWillBeDeleted();
}

//void UbTableModel::objectChanged(UbObservable* changedObject)
//{
//	this->notifyObserversObjectChanged();	
//}
//
//void UbTableModel::objectWillBeDeleted(UbObservable* objectForDeletion)
//{
//	objectForDeletion->removeObserver(this);
//}
