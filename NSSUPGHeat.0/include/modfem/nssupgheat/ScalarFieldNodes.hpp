#ifndef SCALARFIELDNODES_HPP
#define SCALARFIELDNODES_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API ScalarFieldNodes:
	public QObject
{
		Q_OBJECT

	public:
		typedef float real;

		Q_PROPERTY(QByteArray values READ values NOTIFY valuesChanged)

		Q_PROPERTY(int count READ count NOTIFY countChanged)

		ScalarFieldNodes(QObject * parent = nullptr);

		QByteArray values() const;

		int count() const;

	signals:
		void valuesChanged();

		void countChanged();

	private:
		struct Members {
			QByteArray values;
			int count;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // SCALARFIELDNODES_HPP
