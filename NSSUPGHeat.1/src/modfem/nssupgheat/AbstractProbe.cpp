#include "../../../include/modfem/nssupgheat/AbstractProbe.hpp"

namespace modfem {
namespace nssupgheat {

AbstractProbe::AbstractProbe(QObject * parent):
	QObject(parent),
	m(new Members)
{
}

QVector3D AbstractProbe::position() const
{
	return m->position;
}

void AbstractProbe::setPosition(const QVector3D position)
{
	if (m->position != position) {
		m->position = position;
		emit positionChanged();
	}
}

QString AbstractProbe::field() const
{
	return m->field;
}

void AbstractProbe::setField(const QString & field)
{
	if (m->field != field) {
		m->field = field;
		emit fieldChanged();
	}
}

}
}
