#ifndef H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_BOUNDARYCONDITIONSDATA_HPP
#define H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_BOUNDARYCONDITIONSDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace heat {

class MODFEM_HEAT_API BoundaryConditionsData:
	public QObject
{
		Q_OBJECT

	public:
		explicit BoundaryConditionsData(QObject * parent = nullptr);

	public slots:
		void init(int meshId);

	signals:

	private:
		struct Members {
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // MESH_HPP
