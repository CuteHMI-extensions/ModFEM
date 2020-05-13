#include "../../../include/modfem/nssupgheat/BoundaryConditionsData.hpp"

//#include <modfem/mmh_intf.h>

namespace modfem {
namespace nssupgheat {

BoundaryConditionsData::BoundaryConditionsData(QObject * parent):
	QObject(parent),
	m(new Members{})
{
}

void BoundaryConditionsData::init(int meshId)
{
}

}
}
