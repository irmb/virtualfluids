/****************************************************************************
** Smtp meta object code from reading C++ file 'smtp.h'
**
** Created: Mo 22. Nov 19:16:09 2004
**      by: The Qt MOC ($Id: $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "./smtp.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *Smtp::className() const
{
    return "Smtp";
}

QMetaObject *Smtp::metaObj = 0;
static QMetaObjectCleanUp cleanUp_Smtp( "Smtp", &Smtp::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString Smtp::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "Smtp", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString Smtp::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "Smtp", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* Smtp::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QObject::staticMetaObject();
    static const QUMethod slot_0 = {"dnsLookupHelper", 0, 0 };
    static const QUMethod slot_1 = {"readyRead", 0, 0 };
    static const QUMethod slot_2 = {"connected", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "dnsLookupHelper()", &slot_0, QMetaData::Private },
	{ "readyRead()", &slot_1, QMetaData::Private },
	{ "connected()", &slot_2, QMetaData::Private }
    };
    static const QUParameter param_signal_0[] = {
	{ 0, &static_QUType_QString, 0, QUParameter::In }
    };
    static const QUMethod signal_0 = {"status", 1, param_signal_0 };
    static const QMetaData signal_tbl[] = {
	{ "status(const QString&)", &signal_0, QMetaData::Public }
    };
    metaObj = QMetaObject::new_metaobject(
	"Smtp", parentObject,
	slot_tbl, 3,
	signal_tbl, 1,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_Smtp.setMetaObject( metaObj );
    return metaObj;
}

void* Smtp::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "Smtp" ) )
	return this;
    return QObject::qt_cast( clname );
}

// SIGNAL status
void Smtp::status( const QString& t0 )
{
    activate_signal( staticMetaObject()->signalOffset() + 0, t0 );
}

bool Smtp::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: dnsLookupHelper(); break;
    case 1: readyRead(); break;
    case 2: connected(); break;
    default:
	return QObject::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool Smtp::qt_emit( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->signalOffset() ) {
    case 0: status((const QString&)static_QUType_QString.get(_o+1)); break;
    default:
	return QObject::qt_emit(_id,_o);
    }
    return TRUE;
}
#ifndef QT_NO_PROPERTIES

bool Smtp::qt_property( int id, int f, QVariant* v)
{
    return QObject::qt_property( id, f, v);
}

bool Smtp::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
