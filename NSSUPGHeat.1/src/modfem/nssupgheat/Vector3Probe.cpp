#include "../../../include/modfem/nssupgheat/Vector3Probe.hpp"

namespace modfem {
namespace nssupgheat {

modfem::nssupgheat::Vector3Probe::Vector3Probe(QObject * parent):
	AbstractProbe(parent),
	m(new Members{{std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()}})
{
}

QVector3D Vector3Probe::value() const
{
	return m->value;
}

void Vector3Probe::setValue(const QVector3D & value)
{
	if (m->value != value) {
		m->value = value;
		emit valueChanged();
	}
}

}
}
