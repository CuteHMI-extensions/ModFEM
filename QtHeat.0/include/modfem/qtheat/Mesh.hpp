#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_MESH_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_MESH_HPP

#include "internal/common.hpp"
#include "FaceData.hpp"

#include <QObject>

namespace modfem {
namespace qtheat {

class MODFEM_QTHEAT_API Mesh:
	public QObject
{
		Q_OBJECT

	public:
		/// @todo switch beteen active/inactive.
//		Q_PROPERTY(bool active READ active WRITE setActive NOTIFY activeChanged)

		Q_PROPERTY(int nodeCount READ nodeCount NOTIFY nodeCountChanged)

		Q_PROPERTY(QByteArray nodeCoords READ nodeCoords NOTIFY nodeCoordsChanged)

		Q_PROPERTY(FaceData * faceData READ faceData NOTIFY faceDataChanged)

		Q_PROPERTY(QVariantMap nodes READ nodes NOTIFY nodesChanged)

		Q_PROPERTY(QVariantMap triangles READ triangles NOTIFY trianglesChanged)

		Q_PROPERTY(QVariantMap quads READ quads NOTIFY quadsChanged)

		Q_PROPERTY(QVariantMap lines READ lines NOTIFY linesChanged)

		explicit Mesh(QObject * parent = nullptr);

		int nodeCount() const;

		QByteArray nodeCoords() const;

		FaceData * faceData() const;

		QVariantMap nodes() const;

		QVariantMap triangles() const;

		QVariantMap quads() const;

		QVariantMap lines() const;

	public slots:
		void init(int meshId);

	signals:
		void nodeCoordsChanged();

		void nodeCountChanged();

		void faceDataChanged();

		void nodesChanged();

		void trianglesChanged();

		void quadsChanged();

		void linesChanged();

	private slots:
		void setNodeData(int meshId);

		void setNodeCount(int count);

		void update(int meshId);

		void count(int meshId);

	private:
		struct Members {
			int nodeCount;
			QByteArray nodeCoords;	///@todo rename to activeNodeVertices and add inactiveNodeVertices array.
			FaceData * faceData;
			QVariantMap nodes;
			QVariantMap triangles;
			QVariantMap quads;
			QVariantMap lines;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // MESH_HPP
