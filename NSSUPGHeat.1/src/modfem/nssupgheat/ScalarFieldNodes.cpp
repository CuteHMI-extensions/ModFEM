#include "../../../include/modfem/nssupgheat/ScalarFieldNodes.hpp"

namespace modfem {
namespace nssupgheat {

ScalarFieldNodes::ScalarFieldNodes(QObject * parent):
	QObject(parent),
	m(new Members)
{
}

QByteArray ScalarFieldNodes::values() const
{
	return m->values;
}

int ScalarFieldNodes::count() const
{
	return m->count;
}

}
}
