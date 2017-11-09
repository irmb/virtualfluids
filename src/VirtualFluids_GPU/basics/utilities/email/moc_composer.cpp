/****************************************************************************
** Composer meta object code from reading C++ file 'composer.h'
**
** Created: Mo 22. Nov 19:37:52 2004
**      by: The Qt MOC ($Id: $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "./composer.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *Composer::className() const
{
    return "Composer";
}

QMetaObject *Composer::metaObj = 0;
static QMetaObjectCleanUp cleanUp_Composer( "Composer", &Composer::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString Composer::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "Composer", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString Composer::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "Composer", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* Composer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QWidget::staticMetaObject();
    static const QUMethod slot_0 = {"sendMessage", 0, 0 };
    static const QUMethod slot_1 = {"enableSend", 0, 0 };
    static const QUParameter param_slot_2[] = {
	{ "text", &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod slot_2 = {"coutText", 1, param_slot_2 };
    static const QUMethod slot_3 = {"reduceEmailCounter", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "sendMessage()", &slot_0, QMetaData::Private },
	{ "enableSend()", &slot_1, QMetaData::Private },
	{ "coutText(const QString&)", &slot_2, QMetaData::Private },
	{ "reduceEmailCounter()", &slot_3, QMetaData::Private }
    };
    metaObj = QMetaObject::new_metaobject(
	"Composer", parentObject,
	slot_tbl, 4,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_Composer.setMetaObject( metaObj );
    return metaObj;
}

void* Composer::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "Composer" ) )
	return this;
    return QWidget::qt_cast( clname );
}

bool Composer::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: sendMessage(); break;
    case 1: enableSend(); break;
    case 2: coutText((const QString&)static_QUType_QString.get(_o+1)); break;
    case 3: reduceEmailCounter(); break;
    default:
	return QWidget::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool Composer::qt_emit( int _id, QUObject* _o )
{
    return QWidget::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool Composer::qt_property( int id, int f, QVariant* v)
{
    return QWidget::qt_property( id, f, v);
}

bool Composer::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
