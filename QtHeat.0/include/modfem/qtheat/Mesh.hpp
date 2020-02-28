#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_MESH_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_MESH_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace qtheat {

class MODFEM_QTHEAT_API Mesh:
	public QObject
{
		Q_OBJECT

	public:
		Q_PROPERTY(int nodeCount READ nodeCount NOTIFY nodeCountChanged)

		Q_PROPERTY(QByteArray nodeData READ nodeData NOTIFY nodeDataChanged)

		explicit Mesh(QObject * parent = nullptr);

		int nodeCount() const;

		QByteArray nodeData() const;

	public slots:
		void init(int meshId);

	signals:
		void nodeDataChanged();

		void nodeCountChanged();

	private slots:
		void setNodeData(int meshId);

		void setNodeCount(int count);

	private:
		struct Members {
			int nodeCount;
			QByteArray nodeData;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // MESH_HPP
