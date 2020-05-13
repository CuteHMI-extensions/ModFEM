#ifndef H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP
#define H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API ElementData:
	public QObject
{
		Q_OBJECT

	public:
		typedef float real;

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

		void selectAll();

		void selectElement(int index);

	signals:
		void countChanged();

		void nodesChanged();

		void trianglesChanged();

		void quadsChanged();

		void linesChanged();

	private:
		// Dimension of coordinate.
		static constexpr std::size_t COORD_DIM = 3;

		// Size of coordinate in bytes.
		static constexpr std::size_t COORD_SIZE = sizeof(real) * COORD_DIM;

		// Dimension of index.
		static constexpr std::size_t INDEX_DIM = 1;

		// Size of index in bytes.
		static constexpr std::size_t INDEX_SIZE = sizeof(int) * INDEX_DIM;

		// Dimension of normal vector.
		static constexpr std::size_t NORMAL_DIM = 3;

		// Size of normal vector in bytes.
		static constexpr std::size_t NORMAL_SIZE = sizeof(real) * NORMAL_DIM;

		// Dimension of boundary condition.
		static constexpr std::size_t BOUNDARY_CONDITION_DIM = 1;

		// Size of boundary condition in bytes.
		static constexpr std::size_t BOUNDARY_CONDITION_SIZE = sizeof(int) * BOUNDARY_CONDITION_DIM;

		// Dimension of temperature.
		static constexpr std::size_t TEMPERATURE_DIM = 1;

		// Size of temperature scalar in bytes.
		static constexpr std::size_t TEMPERATURE_SIZE = sizeof(real) * TEMPERATURE_DIM;

		void clearArrays();

		void reserveArrays();

		void updateArrays(int meshId);

		void updateArrays(int meshId, int elementId);

		void updateTriangleArrays(int meshId, int elementId);

		void updateProperties();

		void countEntities(int meshId);

		template <std::size_t DIM, typename T>
		void appendRealVector(T * vector,  QByteArray & array);

		template <std::size_t DIM, typename T>
		void appendScalar(T scalar,  QByteArray & array);

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
			QByteArray triangleBoundaryConditions;
			QByteArray quadCoords;
			QByteArray quadIndices;
			QByteArray quadNormals;
			QByteArray quadBoundaryConditions;
			QByteArray quadTriangleCoords;
			QByteArray quadTriangleIndices;
			QByteArray quadTriangleNormals;
			QByteArray quadTriangleBoundaryConditions;
			QByteArray triangleTemperatures;
		};

		cutehmi::MPtr<Members> m;
};

template<std::size_t DIM, typename T>
void ElementData::appendRealVector(T * vector, QByteArray & array)
{
	real realVector[DIM];
	for (std::size_t i = 0; i < DIM; i++)
		realVector[i] = vector[i];

	array.append(reinterpret_cast<char *>(realVector), sizeof(real) * DIM);
}

template<std::size_t DIM, typename T>
void ElementData::appendScalar(T scalar, QByteArray & array)
{
	array.append(reinterpret_cast<char *>(& scalar), sizeof(T) * DIM);
}

}
}

#endif // MESH_HPP
