#include "../../../include/modfem/nssupgheat/ScalarProbe.hpp"

namespace modfem {
namespace nssupgheat {

ScalarProbe::ScalarProbe(QObject * parent):
	AbstractProbe(parent),
	m(new Members{std::numeric_limits<double>::quiet_NaN()})
{
}

double ScalarProbe::value() const
{
	return m->value;
}

void ScalarProbe::setValue(double value)
{
	if (m->value != value) {
		m->value = value;
		emit valueChanged();
	}
}

}
}
