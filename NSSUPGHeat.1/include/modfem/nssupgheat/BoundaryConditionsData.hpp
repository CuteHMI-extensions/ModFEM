#ifndef H_EXTENSIONS_MODFEM_NSSUPGHEAT_0_INCLUDE_MODFEM_NSSUPGHEAT_BOUNDARYCONDITIONSDATA_HPP
#define H_EXTENSIONS_MODFEM_NSSUPGHEAT_0_INCLUDE_MODFEM_NSSUPGHEAT_BOUNDARYCONDITIONSDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API BoundaryConditionsData:
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
