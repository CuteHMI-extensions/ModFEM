#include "../../../include/modfem/nssupgheat/ElementData.hpp"

#include <modfem/mmh_intf.h>
#include <modfem/aph_intf.h>
#include <modfem/pd_ns_supg/pdh_ns_supg_problem.h>
#include <modfem/pd_heat/pdh_heat_problem.h>
#include <modfem/pd_ns_supg_heat/pdh_ns_supg_heat.h>

namespace modfem {
namespace nssupgheat {

ElementData::ElementData(QObject * parent):
	QObject(parent),
	m(new Members{
	0,
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{},
	{}})
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

void ElementData::init(int meshId)
{
	m->meshId = meshId;
	countEntities(meshId);
	reserveArrays();
	selectAll();
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
	} else
		CUTEHMI_WARNING("Attempting to select elements on uninitialized mesh.");
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
//	const int temperatureId = 0;

	int solutionCount = apr_get_nreq(fieldId);
//	CUTEHMI_DEBUG("solutionCount: " << solutionCount);

//	CUTEHMI_DEBUG("apr_get_el_dofs: " << apr_get_el_dofs(fieldId, elementId, 1, elDofs));
	apr_get_el_dofs(fieldId, elementId, 1, elDofs);
	int nodeCount = mmr_el_node_coor(m->meshId, elementId, elNodes, NULL);
	int dofCtr = 0;
	for (int nodeIt = 1; nodeIt <= nodeCount; nodeIt++) {
//		int nrdofs = apr_get_ent_nrdofs(fieldId, APC_VERTEX, nodeIt);
		for (int dofIt = 0; dofIt < solutionCount; dofIt++) {
//			CUTEHMI_DEBUG("nodeIt: " << nodeIt << " nrdofs: " << nrdofs << " elNodes[nodeIt]: " << elNodes[nodeIt] << " elDofs[dofCtr]: " << elDofs[dofCtr++]);
//			CUTEHMI_DEBUG("Insert into " << elNodes[nodeIt] << " value " << elDofs[dofCtr] << " reserved size " << m->nodeTemperatures.size());
			m->nodeTemperatures[elNodes[nodeIt]] = elDofs[dofCtr++];
//			int nodeId = elNodes[nodeIt];
////			mf_check(dof_counter < (sizeof el_dofs) / sizeof(double), "Dof counter (%d) exceeds el_dofs size (%d)",	dof_counter, (sizeof el_dofs) / sizeof(double) );
//			solInfos[node_id].dofs[dofIndex] = elDofs[dofCounter];
//			//Fk-for testing material
//			//solInfos[node_id].dofs[idof]=mmr_el_groupID(mesh_id, el_id);

//			dofCounter++;
		}
	}

	fieldId = pdv_ns_supg_problem.ctrl.field_id;
	solutionCount = apr_get_nreq(fieldId);
//	CUTEHMI_DEBUG("solutionCount: " << solutionCount);
//	CUTEHMI_DEBUG("apr_get_el_dofs: " << apr_get_el_dofs(fieldId, elementId, 1, elDofs));
	apr_get_el_dofs(fieldId, elementId, 1, elDofs);
	nodeCount = mmr_el_node_coor(m->meshId, elementId, elNodes, NULL);
	dofCtr = 0;
	for (int nodeIt = 1; nodeIt <= nodeCount; nodeIt++) {
//		int nrdofs = apr_get_ent_nrdofs(fieldId, APC_VERTEX, nodeIt);
		for (int dofIt = 0; dofIt < solutionCount; dofIt++) {
//			CUTEHMI_DEBUG("nodeIt: " << nodeIt << " nrdofs: " << nrdofs << " elNodes[nodeIt]: " << elNodes[nodeIt] << " elDofs[dofCtr]: " << elDofs[dofCtr++]);
			m->nodeVelocities[elNodes[nodeIt]][0] = elDofs[dofCtr++];
			m->nodeVelocities[elNodes[nodeIt]][1] = elDofs[dofCtr++];
			m->nodeVelocities[elNodes[nodeIt]][2] = elDofs[dofCtr++];
			m->nodePressures[elNodes[nodeIt]] = elDofs[dofCtr++];
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

	for (int i = 1; i <= 3; i++)
		appendScalar(m->nodeTemperatures[indices[i]], m->triangleTemperatures);
}

void ElementData::assignQuadFields(int meshId, int faceId)
{
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	mmr_fa_node_coor(meshId, faceId, indices, nullptr);

	for (int i = 1; i <= 4; i++)
		appendScalar(m->nodeTemperatures[indices[i]], m->quadTemperatures);
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
}

void ElementData::updateFieldProperties()
{
	m->triangleFields["temperatures"] = m->triangleTemperatures;
	emit triangleFieldsChanged();

	m->quadFields["temperatures"] = m->quadTemperatures;
	m->quadFields["triangleTemperatures"] = m->quadTriangleTemperatures;
	emit quadFieldsChanged();
}

void ElementData::clearFieldArrays()
{
	m->triangleTemperatures.clear();
	m->quadTemperatures.clear();
	m->quadTriangleTemperatures.clear();
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
