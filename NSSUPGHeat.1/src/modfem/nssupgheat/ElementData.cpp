#include "../../../include/modfem/nssupgheat/ElementData.hpp"

#include <modfem/mmh_intf.h>
#include <modfem/aph_intf.h>
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>

#include <cmath>

namespace modfem {
namespace nssupgheat {

ElementData::ElementData(QObject * parent):
	QObject(parent),
	m(new Members(this))
{
	m->count["elements"] = 0;
//	m->count["faces"] = 0;

	m->triangles["count"] = 0;
	m->quads["count"] = 0;
	m->quads["triangleCount"] = 0;
	m->lines["count"] = 0;
	m->nodes["count"] = 0;
}

QVariantMap ElementData::count() const
{
	return m->count;
}

QVariantMap ElementData::nodes() const
{
	return m->nodes;
}

QVariantMap ElementData::triangles() const
{
	return m->triangles;
}

QVariantMap ElementData::quads() const
{
	return m->quads;
}

QVariantMap ElementData::lines() const
{
	return m->lines;
}

QVariantMap ElementData::triangleFields() const
{
	return m->triangleFields;
}

QVariantMap ElementData::quadFields() const
{
	return m->quadFields;
}

double ElementData::minTemperature() const
{
	return m->minTemperature;
}

double ElementData::maxTemperature() const
{
	return m->maxTemperature;
}

double ElementData::minPressure() const
{
	return m->minPressure;
}

double ElementData::maxPressure() const
{
	return m->maxPressure;
}

double ElementData::minVelocityMagnitude() const
{
	return m->minVelocityMagnitude;
}

double ElementData::maxVelocityMagnitude() const
{
	return m->maxVelocityMagnitude;
}

QQmlListProperty<AbstractProbe> ElementData::probeList()
{
	return m->probeList;
}

void ElementData::init(int meshId)
{
	m->meshId = meshId;
	clearRecords();
	countEntities(meshId);
	reserveArrays();
	selectAll();
	clearProbes();
	attachProbes();
}

void ElementData::selectAll()
{
	if (m->meshId != 0) {
		clearArrays();
		updateArrays(m->meshId);
		updateProperties();

		clearFieldArrays();
		assignFieldValues();
		updateFieldProperties();
	}
}

void ElementData::selectElement(int elementId)
{
	if (m->meshId != 0) {
		clearArrays();
		updateArrays(m->meshId, elementId);
		updateProperties();

		clearFieldArrays();
		assignFieldValues(elementId);
		updateFieldProperties();
	} else
		CUTEHMI_WARNING("Attempting to select the element on uninitialized mesh.");
}

void ElementData::updateFields()
{
	int elementId = mmr_get_next_act_elem(m->meshId, 0);

	while (elementId != 0) {
		updateFields(elementId);
		elementId = mmr_get_next_act_elem(m->meshId, elementId);
	}

	clearFieldArrays();
	assignFieldValues();
	updateFieldProperties();
}

void ElementData::updateFields(int elementId)
{
	// Read fields.

	double elDofs[APC_DOUBLE_MAXELSD] = {0};
	int elNodes[MMC_MAXELVNO + 1] = {0};

	int fieldId = pdv_heat_problem.ctrl.field_id;
	int solutionCount = apr_get_nreq(fieldId);
	apr_get_el_dofs(fieldId, elementId, 1, elDofs);
	int nodeCount = mmr_el_node_coor(m->meshId, elementId, elNodes, NULL);
	int dofCtr = 0;
	for (int nodeIt = 1; nodeIt <= nodeCount; nodeIt++) {
		for (int dofIt = 0; dofIt < solutionCount; dofIt++) {
			m->nodeTemperatures[elNodes[nodeIt]] = elDofs[dofCtr++];
			maybeSetMinTemperature(static_cast<double>(m->nodeTemperatures[elNodes[nodeIt]]));
			maybeSetMaxTemperature(static_cast<double>(m->nodeTemperatures[elNodes[nodeIt]]));
		}
	}

	fieldId = pdv_ns_supg_problem.ctrl.field_id;
	solutionCount = apr_get_nreq(fieldId);
	apr_get_el_dofs(fieldId, elementId, 1, elDofs);
	nodeCount = mmr_el_node_coor(m->meshId, elementId, elNodes, NULL);
	dofCtr = 0;
	for (int nodeIt = 1; nodeIt <= nodeCount; nodeIt++) {
		for (int dofIt = 0; dofIt < solutionCount; dofIt++) {
			m->nodeVelocities[elNodes[nodeIt]][0] = elDofs[dofCtr++];
			m->nodeVelocities[elNodes[nodeIt]][1] = elDofs[dofCtr++];
			m->nodeVelocities[elNodes[nodeIt]][2] = elDofs[dofCtr++];
			maybeSetMinVelocityMagnitude(norm(m->nodeVelocities[elNodes[nodeIt]]));
			maybeSetMaxVelocityMagnitude(norm(m->nodeVelocities[elNodes[nodeIt]]));
			m->nodePressures[elNodes[nodeIt]] = elDofs[dofCtr++];
			maybeSetMinPressure(static_cast<double>(m->nodePressures[elNodes[nodeIt]]));
			maybeSetMaxPressure(static_cast<double>(m->nodePressures[elNodes[nodeIt]]));
		}
	}
}

void ElementData::assignFieldValues()
{
	int elementId = mmr_get_next_act_elem(m->meshId, 0);
	while (elementId != 0) {
		assignFieldValues(elementId);
		elementId = mmr_get_next_act_elem(m->meshId, elementId);
	}
}

void ElementData::assignFieldValues(int elementId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces in 1st element).

	int faces[FACES_MAX] = {0};
	int orientation[FACES_MAX];

	mmr_el_faces(m->meshId, elementId, faces, orientation);
	CUTEHMI_ASSERT(faces[0] < FACES_MAX, "not enough space reserved for 'faces' array");

	int faceCount = faces[0];
	for (int i = 1; i <= faceCount; i++) {
		switch (mmr_fa_type(m->meshId, faces[i])) {
			case MMC_TRIA:
				assignTriangleFields(m->meshId, faces[i]);
				break;
			case MMC_QUAD:
				assignQuadFields(m->meshId, faces[i]);
				assignQuadTriangleFields(m->meshId, faces[i]);	// Tessellation of quads due to lack of GL_QUADS in modern OpenGL.
				break;
			default:
				CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(m->meshId, faces[i]) << "'");
		}
	}
}

void ElementData::updateProbes()
{
	for (auto it = m->probeScalarFieldNodeMap.begin(); it != m->probeScalarFieldNodeMap.end(); ++it)
		it.key()->setValue(it.value().first->at(it.value().second));
	for (auto it = m->probeVector3FieldNodeMap.begin(); it != m->probeVector3FieldNodeMap.end(); ++it) {
		Vector3FieldContainer * container = it.value().first;
		std::array<freal, 3> value = container->at(it.value().second);
		it.key()->setValue(QVector3D(static_cast<float>(value[0]), static_cast<float>(value[1]), static_cast<float>(value[2])));
	}
}

int ElementData::ProbeListCount(QQmlListProperty<AbstractProbe> * property)
{
	return static_cast<ProbeContainer *>(property->data)->count();
}

AbstractProbe * ElementData::ProbeListAt(QQmlListProperty<AbstractProbe> * property, int index)
{
	return static_cast<ProbeContainer *>(property->data)->value(index);
}

void ElementData::ProbeListClear(QQmlListProperty<AbstractProbe> * property)
{
	ElementData * elementData = static_cast<ElementData *>(property->object);
	for (ProbeContainer::const_iterator it = elementData->probes().begin(); it != elementData->probes().end(); ++it) {
		(*it)->disconnect(elementData);
	}
	static_cast<ProbeContainer *>(property->data)->clear();
}

void ElementData::ProbeListAppend(QQmlListProperty<AbstractProbe> * property, AbstractProbe * value)
{
	static_cast<ProbeContainer *>(property->data)->append(value);
}

const ElementData::ProbeContainer & ElementData::probes() const
{
	return m->probes;
}

void ElementData::clearRecords()
{
	m->minTemperature = std::numeric_limits<double>::max();
	m->maxTemperature = std::numeric_limits<double>::min();
	m->minPressure = std::numeric_limits<double>::max();
	m->maxPressure = std::numeric_limits<double>::min();
	m->minVelocityMagnitude = std::numeric_limits<double>::max();
	m->maxVelocityMagnitude = std::numeric_limits<double>::min();
}

void ElementData::clearProbes()
{
	for (auto probe : m->probes)
		probe->disconnect(this);
	m->probeScalarFieldNodeMap.clear();
	m->probeVector3FieldNodeMap.clear();
}

void ElementData::maybeSetMinTemperature(double temperature)
{
	if (m->minTemperature > temperature) {
		m->minTemperature = temperature;
		emit minTemperatureChanged();
	}
}

void ElementData::maybeSetMaxTemperature(double temperature)
{
	if (m->maxTemperature < temperature) {
		m->maxTemperature = temperature;
		emit maxTemperatureChanged();
	}
}

void ElementData::maybeSetMinPressure(double pressure)
{
	if (m->minPressure > pressure) {
		m->minPressure = pressure;
		emit minPressureChanged();
	}
}

void ElementData::maybeSetMaxPressure(double pressure)
{
	if (m->maxPressure < pressure) {
		m->maxPressure = pressure;
		emit maxPressureChanged();
	}
}

void ElementData::maybeSetMinVelocityMagnitude(double velocityMagnitude)
{
	if (m->minVelocityMagnitude > velocityMagnitude) {
		m->minVelocityMagnitude = velocityMagnitude;
		emit minVelocityMagnitudeChanged();
	}
}

void ElementData::maybeSetMaxVelocityMagnitude(double velocityMagnitude)
{
	if (m->maxVelocityMagnitude < velocityMagnitude) {
		m->maxVelocityMagnitude = velocityMagnitude;
		emit maxVelocityMagnitudeChanged();
	}
}

void ElementData::clearArrays()
{
	m->nodeCoords.clear();
	m->lineCoords.clear();
	m->triangleCoords.clear();
	m->triangleIndices.clear();
	m->triangleNormals.clear();
	m->triangleBoundaryConditions.clear();
	m->quadCoords.clear();
	m->quadIndices.clear();
	m->quadNormals.clear();
	m->quadBoundaryConditions.clear();
	m->quadTriangleCoords.clear();
	m->quadTriangleIndices.clear();
	m->quadTriangleNormals.clear();
	m->quadTriangleBoundaryConditions.clear();
}

void ElementData::reserveArrays()
{
	m->lineCoords.reserve(m->count["lines"].toInt() * 2 * COORD_SIZE);
	m->triangleCoords.reserve(m->count["triangles"].toInt() * 3 * COORD_SIZE);
	m->triangleIndices.reserve(m->count["triangles"].toInt() * 3 * INDEX_SIZE);
	m->triangleNormals.reserve(m->count["triangles"].toInt() * 3 * NORMAL_SIZE);
	m->triangleBoundaryConditions.reserve(m->count["triangles"].toInt() * 3 * BOUNDARY_CONDITION_SIZE);
	m->quadCoords.reserve(m->count["quads"].toInt() * 4 * COORD_SIZE);
	m->quadIndices.reserve(m->count["quads"].toInt() * 4 * INDEX_SIZE);
	m->quadNormals.reserve(m->count["quads"].toInt() * 4 * NORMAL_SIZE);
	m->quadBoundaryConditions.reserve(m->count["quads"].toInt() * 4 * BOUNDARY_CONDITION_SIZE);
	m->quadTriangleCoords.reserve(m->count["quads"].toInt() * 6 * COORD_SIZE);
	m->quadTriangleIndices.reserve(m->count["quads"].toInt() * 6 * INDEX_SIZE);
	m->quadTriangleNormals.reserve(m->count["quads"].toInt() * 6 * NORMAL_SIZE);
	m->quadTriangleBoundaryConditions.reserve(m->count["quads"].toInt() * 6 * BOUNDARY_CONDITION_SIZE);
	m->nodeCoords.reserve(m->count["nodes"].toInt() * COORD_SIZE);

	m->nodeTemperatures.resize(m->count["nodes"].toInt() + 1);
	m->nodeVelocities.resize(m->count["nodes"].toInt() + 1);
	m->nodePressures.resize(m->count["nodes"].toInt() + 1);

	m->triangleTemperatures.reserve(m->count["triangles"].toInt() * 3 * FREAL_SIZE);
	m->quadTemperatures.reserve(m->count["quads"].toInt() * 4 * FREAL_SIZE);
	m->quadTriangleTemperatures.reserve(m->count["quads"].toInt() * 6 * FREAL_SIZE);

	m->trianglePressures.reserve(m->count["triangles"].toInt() * 3 * FREAL_SIZE);
	m->quadPressures.reserve(m->count["quads"].toInt() * 4 * FREAL_SIZE);
	m->quadTrianglePressures.reserve(m->count["quads"].toInt() * 6 * FREAL_SIZE);

	m->triangleVelocities.reserve(m->count["triangles"].toInt() * 3 * FREAL_SIZE);
	m->quadVelocities.reserve(m->count["quads"].toInt() * 4 * FREAL_SIZE);
	m->quadTriangleVelocities.reserve(m->count["quads"].toInt() * 6 * FREAL_SIZE);
}

void ElementData::updateArrays(int meshId)
{
	int elementIndex = mmr_get_next_act_elem(meshId, 0);

	while (elementIndex != 0) {
		updateArrays(meshId, elementIndex);
		elementIndex = mmr_get_next_act_elem(meshId, elementIndex);
	}

	double coords[3];
	for (int i = 1; i <= mmr_get_max_node_id(meshId); i++) {
		if (mmr_node_status(meshId, i) == MMC_ACTIVE) {
			mmr_node_coor(meshId, i, coords);
			appendAsRealVector<COORD_DIM>(coords, m->nodeCoords);
		}
	}
}

void ElementData::updateArrays(int meshId, int elementId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces in 1st element).

	int faces[FACES_MAX] = {0};
	int orientation[FACES_MAX];

	mmr_el_faces(meshId, elementId, faces, orientation);
	CUTEHMI_ASSERT(faces[0] < FACES_MAX, "not enough space reserved for 'faces' array");

	int faceCount = faces[0];
	for (int i = 1; i <= faceCount; i++) {
		switch (mmr_fa_type(meshId, faces[i])) {
			case MMC_TRIA:
				updateTriangleArrays(meshId, faces[i]);
				break;
			case MMC_QUAD:
				updateQuadArrays(meshId, faces[i]);
				updateQuadTriangleArrays(meshId, faces[i]);	// Tessellation of quads due to lack of GL_QUADS in modern OpenGL.
				break;
			default:
				CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faces[i]) << "'");
		}
	}
}

void ElementData::updateTriangleArrays(int meshId, int faceId)
{
	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
	double vector3[3];	// Space required for x, y, z coordinates.
	//</ModFEM.Mesh-1.workaround>
	int indices[4];		// Triangle faces require 4 indices (first element contains amount of coordinates).

	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// mmr_fa_node_coor(meshId, faceId, indices, coords);
	// triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 3);
	mmr_fa_node_coor(meshId, faceId, indices, nullptr);

	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	mmr_node_coor(meshId, indices[2], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[3], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	//</ModFEM.Mesh-1.workaround>

	appendAsRealVector<INDEX_DIM * 3>(indices + 1, m->triangleIndices);	// First element contains amount of coordinates.

	// Append normal to each vertex.
	double area;		// Face area.
	mmr_fa_area(meshId, faceId, & area, vector3);
	for (int i = 0; i < 3; i++)
		appendAsRealVector<NORMAL_DIM>(vector3, m->triangleNormals);

	// Append boundary condition flags to each vertex.
	int bc = mmr_fa_bc(meshId, faceId);
	for (int i = 0; i < 3; i++)
		appendScalar(bc, m->triangleBoundaryConditions);
}

void ElementData::updateQuadArrays(int meshId, int faceIndex)
{
	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
	double vector3[3];	// Space required for x, y, z coordinates.
	//</ModFEM.Mesh-1.workaround>
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).

	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// mmr_fa_node_coor(meshId, faces[i], indices, coords);
	// quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 4);
	mmr_fa_node_coor(meshId, faceIndex, indices, nullptr);

	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	mmr_node_coor(meshId, indices[2], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[3], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[4], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadCoords);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->lineCoords);
	//</ModFEM.Mesh-1.workaround>

	appendAsRealVector<INDEX_DIM * 4>(indices + 1, m->quadIndices);	// First element contains amount of coordinates.

	// Append normal to each vertex.
	double area;		// Face area.
	mmr_fa_area(meshId, faceIndex, & area, vector3);
	for (int i = 0; i < 4; i++)
		appendAsRealVector<NORMAL_DIM>(vector3, m->quadNormals);

	// Append boundary condition flags to each vertex.
	int bc = mmr_fa_bc(meshId, faceIndex);
	for (int i = 0; i < 4; i++)
		appendScalar(bc, m->quadBoundaryConditions);
}

void ElementData::updateQuadTriangleArrays(int meshId, int faceIndex)
{
	// Tessellation of quads due to lack of GL_QUADS in modern OpenGL.

	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
	double vector3[3];	// Space required for x, y, z coordinates.
	//</ModFEM.Mesh-1.workaround>
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).

	mmr_fa_node_coor(meshId, faceIndex, indices, nullptr);

	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 1, m->quadTriangleIndices);

	mmr_node_coor(meshId, indices[2], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 2, m->quadTriangleIndices);

	mmr_node_coor(meshId, indices[3], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 3, m->quadTriangleIndices);

	mmr_node_coor(meshId, indices[4], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 4, m->quadTriangleIndices);

	mmr_node_coor(meshId, indices[3], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 3, m->quadTriangleIndices);

	mmr_node_coor(meshId, indices[1], vector3);
	appendAsRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
	appendAsRealVector<INDEX_DIM>(indices + 1, m->quadTriangleIndices);

	// Append normal to each vertex.
	double area;		// Face area.
	mmr_fa_area(meshId, faceIndex, & area, vector3);
	for (int i = 0; i < 6; i++)
		appendAsRealVector<NORMAL_DIM>(vector3, m->quadTriangleNormals);

	// Append boundary condition flags to each vertex.
	int bc = mmr_fa_bc(meshId, faceIndex);
	for (int i = 0; i < 6; i++)
		appendScalar(bc, m->quadTriangleBoundaryConditions);
}

void ElementData::assignTriangleFields(int meshId, int faceId)
{
	int indices[4];		// Triangle faces require 4 indices (first element contains amount of coordinates).
	mmr_fa_node_coor(meshId, faceId, indices, nullptr);

	for (int i = 1; i <= 3; i++) {
		appendScalar(m->nodeTemperatures[indices[i]], m->triangleTemperatures);
		appendScalar(m->nodePressures[indices[i]], m->trianglePressures);
		appendScalar(norm(m->nodeVelocities[indices[i]]), m->triangleVelocities);
	}
}

void ElementData::assignQuadFields(int meshId, int faceId)
{
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	mmr_fa_node_coor(meshId, faceId, indices, nullptr);

	for (int i = 1; i <= 4; i++) {
		appendScalar(m->nodeTemperatures[indices[i]], m->quadTemperatures);
		appendScalar(m->nodePressures[indices[i]], m->quadPressures);
		appendScalar(norm(m->nodeVelocities[indices[i]]), m->quadVelocities);
	}
}

void ElementData::assignQuadTriangleFields(int meshId, int faceId)
{
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	mmr_fa_node_coor(meshId, faceId, indices, nullptr);

	appendScalar(m->nodeTemperatures[indices[1]], m->quadTriangleTemperatures);
	appendScalar(m->nodeTemperatures[indices[2]], m->quadTriangleTemperatures);
	appendScalar(m->nodeTemperatures[indices[3]], m->quadTriangleTemperatures);
	appendScalar(m->nodeTemperatures[indices[4]], m->quadTriangleTemperatures);
	appendScalar(m->nodeTemperatures[indices[3]], m->quadTriangleTemperatures);
	appendScalar(m->nodeTemperatures[indices[1]], m->quadTriangleTemperatures);

	appendScalar(m->nodePressures[indices[1]], m->quadTrianglePressures);
	appendScalar(m->nodePressures[indices[2]], m->quadTrianglePressures);
	appendScalar(m->nodePressures[indices[3]], m->quadTrianglePressures);
	appendScalar(m->nodePressures[indices[4]], m->quadTrianglePressures);
	appendScalar(m->nodePressures[indices[3]], m->quadTrianglePressures);
	appendScalar(m->nodePressures[indices[1]], m->quadTrianglePressures);

	appendScalar(norm(m->nodeVelocities[indices[1]]), m->quadTriangleVelocities);
	appendScalar(norm(m->nodeVelocities[indices[2]]), m->quadTriangleVelocities);
	appendScalar(norm(m->nodeVelocities[indices[3]]), m->quadTriangleVelocities);
	appendScalar(norm(m->nodeVelocities[indices[4]]), m->quadTriangleVelocities);
	appendScalar(norm(m->nodeVelocities[indices[3]]), m->quadTriangleVelocities);
	appendScalar(norm(m->nodeVelocities[indices[1]]), m->quadTriangleVelocities);
}

void ElementData::updateFieldProperties()
{
	m->triangleFields["temperatures"] = m->triangleTemperatures;
	m->triangleFields["pressures"] = m->trianglePressures;
	m->triangleFields["velocities"] = m->triangleVelocities;
	emit triangleFieldsChanged();

	m->quadFields["temperatures"] = m->quadTemperatures;
	m->quadFields["triangleTemperatures"] = m->quadTriangleTemperatures;
	m->quadFields["pressures"] = m->quadPressures;
	m->quadFields["trianglePressures"] = m->quadTrianglePressures;
	m->quadFields["velocities"] = m->quadVelocities;
	m->quadFields["triangleVelocities"] = m->quadTriangleVelocities;
	emit quadFieldsChanged();
}

void ElementData::clearFieldArrays()
{
	m->triangleTemperatures.clear();
	m->quadTemperatures.clear();
	m->quadTriangleTemperatures.clear();
	m->trianglePressures.clear();
	m->quadPressures.clear();
	m->quadTrianglePressures.clear();
	m->triangleVelocities.clear();
	m->quadVelocities.clear();
	m->quadTriangleVelocities.clear();
}

void ElementData::updateProperties()
{
	m->triangles["coords"] = m->triangleCoords;
	m->triangles["count"] = m->triangleCoords.count() / (static_cast<int>(COORD_SIZE) * 3);
	m->triangles["normals"] = m->triangleNormals;
	m->triangles["boundaryConditions"] = m->triangleBoundaryConditions;
	m->quads["coords"] = m->quadCoords;
	m->quads["count"] = m->quadCoords.count() / (static_cast<int>(COORD_SIZE) * 4);
	m->quads["normals"] = m->quadNormals;
	m->quads["boundaryConditions"] = m->quadBoundaryConditions;
	m->quads["triangleCoords"] = m->quadTriangleCoords;
	m->quads["triangleCount"] = m->quadTriangleCoords.count() / (static_cast<int>(COORD_SIZE) * 3);
	m->quads["triangleNormals"] = m->quadTriangleNormals;
	m->quads["triangleBoundaryConditions"] = m->quadTriangleBoundaryConditions;
	m->lines["coords"] = m->lineCoords;
	m->lines["count"] = m->lineCoords.count() / (static_cast<int>(COORD_SIZE) * 2);
	m->nodes["coords"] = m->nodeCoords;
	m->nodes["count"] = m->nodeCoords.count() / static_cast<int>(COORD_SIZE);

	emit trianglesChanged();
	emit quadsChanged();
	emit linesChanged();
	emit nodesChanged();
}

void ElementData::countEntities(int meshId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces as 1st element).

	int tetraCount = 0;
	int prismCount = 0;
	int brickCount = 0;
	int triangleCount = 0;
	int quadCount = 0;
	int otherCount = 0;
	int lineCount = 0;

	int elementId = mmr_get_next_act_elem(meshId, 0);
	while (elementId != 0) {
		int faces[FACES_MAX] = {0};
		int orientation[FACES_MAX];
		switch (mmr_el_type(meshId, elementId)) {
			case MMC_TETRA:
				tetraCount++;
				break;
			case MMC_PRISM:
				prismCount++;
				break;
			case MMC_BRICK:
				brickCount++;
				break;
			default:
				CUTEHMI_CRITICAL("Unsupported element type: '" << mmr_el_type(meshId, elementId) << "'");
		}
		mmr_el_faces(meshId, elementId, faces, orientation);
		CUTEHMI_ASSERT(faces[0] < FACES_MAX, "not enough space reserved for 'faces' array");

		int faceCount = faces[0];
		for (int i = 1; i <= faceCount; i++) {
			switch (mmr_fa_type(meshId, faces[i])) {
				case MMC_TRIA:
					triangleCount++;
					lineCount += 3;
					break;
				case MMC_QUAD:
					quadCount++;
					lineCount += 4;
					break;
				default:
					otherCount++;
					CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faces[i]) << "'");
			}
		}

		elementId = mmr_get_next_act_elem(meshId, elementId);
	}

	m->count["elements"] = mmr_get_max_elem_id(meshId);
	m->count["tetrahedrons"] = tetraCount;
	m->count["prisms"] = prismCount;
	m->count["bricks"] = brickCount;
	m->count["triangles"] = triangleCount;
	m->count["quads"] = quadCount;
	m->count["faces"] = otherCount + triangleCount + quadCount;
	m->count["edges"] = mmr_get_nr_edge(meshId);
	m->count["nodes"] = mmr_get_nr_node(meshId);
	m->count["lines"] = lineCount;
	emit countChanged();
}

void ElementData::attachProbes()
{
	for (auto probe : m->probes) {
		// Find closest node for each of the probes (in futre approximation may be used).
		findClosestNode(probe);
		connect(probe, & AbstractProbe::positionChanged, this, [this, probe]() {
			findClosestNode(probe);
		});
	}
}

void ElementData::findClosestNode(AbstractProbe * probe)
{
	static constexpr float INITIAL_DOT_PRODUCT = std::numeric_limits<float>::max();

	int elNodes[MMC_MAXELVNO + 1] = {0};
	double coords[3];
	float dotProduct = INITIAL_DOT_PRODUCT;
	int nodeId = 0;
	int elementId = mmr_get_next_act_elem(m->meshId, 0);
	while (elementId != 0) {
		int nodeCount = mmr_el_node_coor(m->meshId, elementId, elNodes, NULL);
		for (int nodeIt = 1; nodeIt <= nodeCount; nodeIt++) {
			// Choose node which has least distance.
			mmr_node_coor(m->meshId, elNodes[nodeIt], coords);
			float newDotProduct = (probe->position() - QVector3D(static_cast<float>(coords[0]), static_cast<float>(coords[1]), static_cast<float>(coords[2]))).lengthSquared(); // Unfortunately QVector3D uses float, but it should be OK - pick does not be to super-precise.
			if (newDotProduct < dotProduct) {
				nodeId = elNodes[nodeIt];
				dotProduct = newDotProduct;
			}
		}
		elementId = mmr_get_next_act_elem(m->meshId, elementId);
	}
	if (dotProduct != INITIAL_DOT_PRODUCT) {
		if (probe->field() == "temperature") {
			ScalarProbe * scalarProbe = qobject_cast<ScalarProbe *>(probe);
			if (scalarProbe)
				m->probeScalarFieldNodeMap.insert(scalarProbe, std::make_pair(& m->nodeTemperatures, nodeId));
			else
				CUTEHMI_CRITICAL("Probe " << probe <<  " type is inadequate for scalar field '" << probe->field() << "'.");
		} else if (probe->field() == "pressure") {
			ScalarProbe * scalarProbe = qobject_cast<ScalarProbe *>(probe);
			if (scalarProbe)
				m->probeScalarFieldNodeMap.insert(scalarProbe, std::make_pair(& m->nodePressures, nodeId));
			else
				CUTEHMI_CRITICAL("Probe " << probe <<  " type is inadequate for scalar field '" << probe->field() << "'.");
		} else if (probe->field() == "velocity") {
			Vector3Probe * vector3Probe = qobject_cast<Vector3Probe *>(probe);
			if (vector3Probe)
				m->probeVector3FieldNodeMap.insert(vector3Probe, std::make_pair(& m->nodeVelocities, nodeId));
			else
				CUTEHMI_CRITICAL("Probe " << probe <<  " type is inadequate for 3D vector field '" << probe->field() << "'.");
		} else
			CUTEHMI_CRITICAL("Unrecognized field name '" << probe->field() << "'.");
	}
}

ElementData::freal ElementData::norm(const std::array<freal, 3> & vector)
{
	return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
}

template<std::size_t DIM, typename T>
void ElementData::appendVector(T * vector, QByteArray & array)
{
	array.append(reinterpret_cast<char *>(vector), sizeof(T) * DIM);
}

template<std::size_t DIM, typename T>
void ElementData::appendAsRealVector(T * vector, QByteArray & array)
{
	greal realVector[DIM];
	for (std::size_t i = 0; i < DIM; i++)
		realVector[i] = vector[i];

	array.append(reinterpret_cast<char *>(realVector), sizeof(greal) * DIM);
}

template<typename T>
void ElementData::appendScalar(T scalar, QByteArray & array)
{
	array.append(reinterpret_cast<char *>(& scalar), sizeof(T));
}

template<std::size_t DIM, typename T>
void ElementData::insertVector(int pos, T * vector, QByteArray & array)
{
	array.insert(pos * sizeof(T) * DIM, reinterpret_cast<char *>(vector), sizeof(T) * DIM);
}

template<std::size_t DIM, typename T>
void ElementData::insertAsRealVector(int pos, T * vector, QByteArray & array)
{
	greal realVector[DIM];
	for (std::size_t i = 0; i < DIM; i++)
		realVector[i] = vector[i];

	array.insert(pos * sizeof(greal) * DIM, reinterpret_cast<char *>(realVector), sizeof(greal) * DIM);
}

template<typename T>
void ElementData::insertScalar(int pos, T scalar, QByteArray & array)
{
	CUTEHMI_DEBUG("insertScalar scalar: " << scalar << " pos: " << pos);	//temp

	array.insert(pos * sizeof(T), reinterpret_cast<char *>(& scalar), sizeof(T));

	//temp
	QByteArray::const_iterator it = array.begin();
	it += pos * sizeof(T);
	T value = *reinterpret_cast<const double *>(it);
	CUTEHMI_DEBUG("recovering Scalar scalar: " << value);
	//endtemp
}

}
}
