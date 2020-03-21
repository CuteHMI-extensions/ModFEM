#ifndef H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP
#define H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace heat {

class MODFEM_HEAT_API ElementData:
	public QObject
{
		Q_OBJECT

	public:
		/// @todo switch beteen active/inactive.

		Q_PROPERTY(QVariantMap count READ count NOTIFY countChanged)

		Q_PROPERTY(QVariantMap nodes READ nodes NOTIFY nodesChanged)

		Q_PROPERTY(QVariantMap triangles READ triangles NOTIFY trianglesChanged)

		// Note: GL_QUADS have been removed from OpenGL, so quads offer 'triangleCoords' and 'triangleNormals' mapping.
		Q_PROPERTY(QVariantMap quads READ quads NOTIFY quadsChanged)

		Q_PROPERTY(QVariantMap lines READ lines NOTIFY linesChanged)

		explicit ElementData(QObject * parent = nullptr);

		QVariantMap count() const;

		QVariantMap nodes() const;

		QVariantMap triangles() const;

		QVariantMap quads() const;

		QVariantMap lines() const;

	public slots:
		void init(int meshId);

//		void selectElement(int index);

	signals:
		void countChanged();

		void nodesChanged();

		void trianglesChanged();

		void quadsChanged();

		void linesChanged();

	private:
		// Size of coordinate in bytes.
		static constexpr int COORD_SIZE = sizeof(double) * 3;

		// Size of index in bytes.
		static constexpr int INDEX_SIZE = sizeof(int);

		// Size of normal vector in bytes.
		static constexpr int NORMAL_SIZE = sizeof(double) * 3;

		void clearArrays();

		void reserveArrays();

		void updateArrays(int meshId);

		void updateArrays(int meshId, int elementId);

		void updateProperties();

		void countEntities(int meshId);

		struct Members {
			int meshId;
			QVariantMap count;
			QVariantMap nodes;
			QVariantMap triangles;
			QVariantMap quads;
			QVariantMap lines;
			QByteArray nodeCoords;
			QByteArray lineCoords;
			QByteArray triangleCoords;
			QByteArray triangleIndices;
			QByteArray triangleNormals;
			QByteArray quadCoords;
			QByteArray quadIndices;
			QByteArray quadNormals;
			QByteArray quadTriangleCoords;
			QByteArray quadTriangleIndices;
			QByteArray quadTriangleNormals;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // MESH_HPP
