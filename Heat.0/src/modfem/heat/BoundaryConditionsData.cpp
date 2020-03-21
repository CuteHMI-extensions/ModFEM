#include "../../../include/modfem/heat/BoundaryConditionsData.hpp"

#include <modfem/mmh_intf.h>

namespace modfem {
namespace heat {

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
