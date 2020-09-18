#ifndef H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP
#define H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_ELEMENTDATA_HPP

#include "internal/common.hpp"
#include "AbstractProbe.hpp"
#include "ScalarProbe.hpp"
#include "Vector3Probe.hpp"

#include <QObject>
#include <QQmlListProperty>
#include <QMultiMap>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API ElementData:
	public QObject
{
		Q_OBJECT

	public:
		typedef float greal;	// Graphics real type.

		Q_PROPERTY(QVariantMap count READ count NOTIFY countChanged)

		Q_PROPERTY(QVariantMap nodes READ nodes NOTIFY nodesChanged)

		Q_PROPERTY(QVariantMap triangles READ triangles NOTIFY trianglesChanged)

		// Note: GL_QUADS have been removed from OpenGL, so quads offer 'triangleCoords' and 'triangleNormals' mapping.
		Q_PROPERTY(QVariantMap quads READ quads NOTIFY quadsChanged)

		Q_PROPERTY(QVariantMap lines READ lines NOTIFY linesChanged)

		Q_PROPERTY(QVariantMap triangleFields READ triangleFields NOTIFY triangleFieldsChanged)

		// Note: GL_QUADS have been removed from OpenGL, so quadFields offer 'triangleTemperatures' mapping.
		Q_PROPERTY(QVariantMap quadFields READ quadFields NOTIFY quadFieldsChanged)

		Q_PROPERTY(QVariantMap capFields READ capFields NOTIFY capFieldsChanged)

		Q_PROPERTY(QVariantMap caps READ caps NOTIFY capsChanged)

		Q_PROPERTY(double minTemperature READ minTemperature NOTIFY minTemperatureChanged)

		Q_PROPERTY(double maxTemperature READ maxTemperature NOTIFY maxTemperatureChanged)

		Q_PROPERTY(double minPressure READ minPressure NOTIFY minPressureChanged)

		Q_PROPERTY(double maxPressure READ maxPressure NOTIFY maxPressureChanged)

		Q_PROPERTY(double minVelocityMagnitude READ minVelocityMagnitude NOTIFY minVelocityMagnitudeChanged)

		Q_PROPERTY(double maxVelocityMagnitude READ maxVelocityMagnitude NOTIFY maxVelocityMagnitudeChanged)

		Q_PROPERTY(QQmlListProperty<modfem::nssupgheat::AbstractProbe> probes READ probeList)

		explicit ElementData(QObject * parent = nullptr);

		QVariantMap count() const;

		QVariantMap nodes() const;

		QVariantMap triangles() const;

		QVariantMap quads() const;

		QVariantMap lines() const;

		QVariantMap triangleFields() const;

		QVariantMap quadFields() const;

		QVariantMap capFields() const;

		QVariantMap caps() const;

		double minTemperature() const;

		double maxTemperature() const;

		double minPressure() const;

		double maxPressure() const;

		double minVelocityMagnitude() const;

		double maxVelocityMagnitude() const;

		QQmlListProperty<AbstractProbe> probeList();

	public slots:
		void init(int meshId);

		void selectAll();

		void selectElement(int elementId);

		void updateFields();

		void updateFields(int elementId);

		void assignFieldValues();

		void assignFieldValues(int elementId);

		void assignCapFields();

		void updateProbes();

		void clip(const QVector4D & plane);

	signals:
		void countChanged();

		void nodesChanged();

		void trianglesChanged();

		void quadsChanged();

		void capsChanged();

		void linesChanged();

		void triangleFieldsChanged();

		void quadFieldsChanged();

		void capFieldsChanged();

		void minTemperatureChanged();

		void maxTemperatureChanged();

		void minPressureChanged();

		void maxPressureChanged();

		void minVelocityMagnitudeChanged();

		void maxVelocityMagnitudeChanged();

	private:
		typedef double freal;	// Field real type.

		typedef QVector<freal> ScalarFieldContainer;
		typedef QVector<std::array<freal, 3>> Vector3FieldContainer;
		typedef QList<AbstractProbe *> ProbeContainer;
		typedef QHash<ScalarProbe *, std::pair<ScalarFieldContainer *, int>> ProbeScalarFieldNodeMap;
		typedef QHash<Vector3Probe *, std::pair<Vector3FieldContainer *, int>> ProbeVector3FieldNodeMap;
		typedef QVector<std::pair<int, QVector3D>> CapIntersectionsContainer;

		// Dimension of coordinate.
		static constexpr std::size_t COORD_DIM = 3;

		// Size of coordinate in bytes.
		static constexpr std::size_t COORD_SIZE = sizeof(greal) * COORD_DIM;

		// Dimension of index.
		static constexpr std::size_t INDEX_DIM = 1;

		// Size of index in bytes.
		static constexpr std::size_t INDEX_SIZE = sizeof(int) * INDEX_DIM;

		// Dimension of normal vector.
		static constexpr std::size_t NORMAL_DIM = 3;

		// Size of normal vector in bytes.
		static constexpr std::size_t NORMAL_SIZE = sizeof(greal) * NORMAL_DIM;

		// Dimension of boundary condition.
		static constexpr std::size_t BOUNDARY_CONDITION_DIM = 1;

		// Size of boundary condition in bytes.
		static constexpr std::size_t BOUNDARY_CONDITION_SIZE = sizeof(int) * BOUNDARY_CONDITION_DIM;

		// Size of field value in bytes.
		static constexpr std::size_t FREAL_SIZE = sizeof(freal);

		static int ProbeListCount(QQmlListProperty<AbstractProbe> * property);

		static AbstractProbe * ProbeListAt(QQmlListProperty<AbstractProbe> * property, int index);

		static void ProbeListClear(QQmlListProperty<AbstractProbe> * property);

		static void ProbeListAppend(QQmlListProperty<AbstractProbe> * property, AbstractProbe * value);

		const ProbeContainer & probes() const;

		void clearRecords();

		void clearProbes();

		void maybeSetMinTemperature(double temperature);

		void maybeSetMaxTemperature(double temperature);

		void maybeSetMinPressure(double pressure);

		void maybeSetMaxPressure(double pressure);

		void maybeSetMinVelocityMagnitude(double velocityMagnitude);

		void maybeSetMaxVelocityMagnitude(double velocityMagnitude);

		void clearArrays();

		void reserveArrays();

		void updateArrays(int meshId);

		void updateArrays(int meshId, int elementId);

		void updateTriangleArrays(int meshId, int faceId);

		void updateQuadArrays(int meshId, int faceIndex);

		void updateQuadTriangleArrays(int meshId, int faceIndex);

		void assignTriangleFields(int meshId, int faceId);

		void assignQuadFields(int meshId, int faceId);

		void assignQuadTriangleFields(int meshId, int faceId);

		void updateFieldProperties();

		void updateCapFieldProperties();

		void clearFieldArrays();

		void clearCapFieldArrays();

		void updateProperties();

		void countEntities(int meshId);

		void attachProbes();

		void findClosestNode(AbstractProbe * probe);

		template <std::size_t DIM, typename T>
		void appendVector(T * vector,  QByteArray & array);

		template <std::size_t DIM, typename T>
		void appendAsRealVector(T * vector,  QByteArray & array);

		void appendAsRealVector(const QVector3D & vector,  QByteArray & array);

		template <typename T>
		void appendScalar(T scalar,  QByteArray & array);

		template <std::size_t DIM, typename T>
		void insertVector(int pos, T * vector,  QByteArray & array);

		template <std::size_t DIM, typename T>
		void insertAsRealVector(int pos, T * vector,  QByteArray & array);

		template <typename T>
		void insertScalar(int pos, T scalar,  QByteArray & array);

		freal norm(const std::array<freal, 3> & vector);

		struct Members {
			int meshId;
			QVariantMap count;
			QVariantMap nodes;
			QVariantMap triangles;
			QVariantMap quads;
			QVariantMap lines;
			QVariantMap triangleFields;
			QVariantMap quadFields;
			QVariantMap caps;
			QVariantMap capFields;
			QByteArray capsCoords;
			QByteArray capsNormals;
			CapIntersectionsContainer capIntersections;
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
			ScalarFieldContainer nodeTemperatures;
			Vector3FieldContainer nodeVelocities;
			ScalarFieldContainer nodePressures;
			QByteArray triangleTemperatures;
			QByteArray quadTemperatures;
			QByteArray quadTriangleTemperatures;
			QByteArray capTemperatures;
			QByteArray trianglePressures;
			QByteArray quadPressures;
			QByteArray quadTrianglePressures;
			QByteArray capPressures;
			QByteArray triangleVelocityMagnitudes;
			QByteArray quadVelocities;
			QByteArray quadTriangleVelocityMagnitudes;
			QByteArray capVelocityMagnitudes;
			double minTemperature;
			double maxTemperature;
			double minPressure;
			double maxPressure;
			double minVelocityMagnitude;
			double maxVelocityMagnitude;
			ProbeContainer probes;
			QQmlListProperty<AbstractProbe> probeList;
			ProbeScalarFieldNodeMap probeScalarFieldNodeMap;
			ProbeVector3FieldNodeMap probeVector3FieldNodeMap;

			Members(ElementData * p_parent):
				meshId(0),
				minTemperature(std::numeric_limits<double>::max()),
				maxTemperature(std::numeric_limits<double>::lowest()),
				minPressure(std::numeric_limits<double>::max()),
				maxPressure(std::numeric_limits<double>::lowest()),
				minVelocityMagnitude(std::numeric_limits<double>::max()),
				maxVelocityMagnitude(std::numeric_limits<double>::lowest()),
				probeList(p_parent, & probes, & ElementData::ProbeListAppend, & ElementData::ProbeListCount, & ElementData::ProbeListAt, & ElementData::ProbeListClear)
			{
			}
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // MESH_HPP
