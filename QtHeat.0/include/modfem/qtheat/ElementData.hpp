#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_ELEMENTDATA_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_ELEMENTDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace qtheat {

class MODFEM_QTHEAT_API ElementData:
	public QObject
{
		Q_OBJECT

	public:
		explicit ElementData(QObject * parent = nullptr);

	signals:
};

}
}

#endif // ELEMENTDATA_HPP
