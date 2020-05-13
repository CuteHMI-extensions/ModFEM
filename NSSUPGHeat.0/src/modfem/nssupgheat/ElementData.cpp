#include "../../../include/modfem/nssupgheat/ElementData.hpp"

#include <modfem/mmh_intf.h>

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
	} else
		CUTEHMI_WARNING("Attempting to select elements on uninitialized mesh.");
}

void ElementData::selectElement(int index)
{
	if (m->meshId != 0) {
		clearArrays();
		updateArrays(m->meshId, index);
		updateProperties();
	} else
		CUTEHMI_WARNING("Attempting to select the element on uninitialized mesh.");
}

void ElementData::clearArrays()
{
	m->nodeCoords.clear();
	m->lineCoords.clear();
	m->triangleCoords.clear();
	m->triangleIndices.clear();
	m->triangleNormals.clear();
	m->triangleBoundaryConditions.clear();
	m->triangleTemperatures.clear();
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
	m->triangleTemperatures.reserve(m->nodes["count"].toInt() * 3 * TEMPERATURE_SIZE);
	m->quadCoords.reserve(m->count["quads"].toInt() * 4 * COORD_SIZE);
	m->quadIndices.reserve(m->count["quads"].toInt() * 4 * INDEX_SIZE);
	m->quadNormals.reserve(m->count["quads"].toInt() * 4 * NORMAL_SIZE);
	m->quadBoundaryConditions.reserve(m->count["quads"].toInt() * 4 * BOUNDARY_CONDITION_SIZE);
	m->quadTriangleCoords.reserve(m->count["quads"].toInt() * 6 * COORD_SIZE);
	m->quadTriangleIndices.reserve(m->count["quads"].toInt() * 6 * INDEX_SIZE);
	m->quadTriangleNormals.reserve(m->count["quads"].toInt() * 6 * NORMAL_SIZE);
	m->quadTriangleBoundaryConditions.reserve(m->count["quads"].toInt() * 6 * BOUNDARY_CONDITION_SIZE);
	m->nodeCoords.reserve(m->nodes["count"].toInt() * COORD_SIZE);
}

void ElementData::updateArrays(int meshId)
{
	int elementId = mmr_get_next_act_elem(meshId, 0);

	while (elementId != 0) {
		updateArrays(meshId, elementId);
		elementId = mmr_get_next_act_elem(meshId, elementId);
	}

	double coords[3];
	for (int i = 1; i <= mmr_get_max_node_id(meshId); i++) {
		if (mmr_node_status(meshId, i) == MMC_ACTIVE) {
			mmr_node_coor(meshId, i, coords);
			appendRealVector<COORD_DIM>(coords, m->nodeCoords);
		}
	}
}

void ElementData::updateArrays(int meshId, int elementId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces in 1st element).

	int faces[FACES_MAX] = {0};
	int orientation[FACES_MAX];
	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
	double vector3[3];	// Space required for x, y, z coordinates.
	//</ModFEM.Mesh-1.workaround>
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	double area;		// Face area.
	int bc;				// Boundary condition flag.

	mmr_el_faces(meshId, elementId, faces, orientation);
	CUTEHMI_ASSERT(faces[0] < FACES_MAX, "not enough space reserved for 'faces' array");

	int faceCount = faces[0];
	for (int i = 1; i <= faceCount; i++) {
		switch (mmr_fa_type(meshId, faces[i])) {
			case MMC_TRIA:
				//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
				// Instead of:
				// mmr_fa_node_coor(meshId, faceId, indices, coords);
				// triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 3);
				mmr_fa_node_coor(meshId, faces[i], indices, nullptr);

				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				mmr_node_coor(meshId, indices[2], vector3);
				appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				mmr_node_coor(meshId, indices[3], vector3);
				appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				//</ModFEM.Mesh-1.workaround>

				appendRealVector<INDEX_DIM * 3>(indices + 1, m->triangleIndices);	// First element contains amount of coordinates.

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, vector3);
				for (int i = 0; i < 3; i++)
					appendRealVector<NORMAL_DIM>(vector3, m->triangleNormals);

				// Append boundary condition flags to each vertex.

				bc = mmr_fa_bc(meshId, faces[i]);
				for (int i = 0; i < 3; i++)
					appendScalar<BOUNDARY_CONDITION_DIM>(bc, m->triangleBoundaryConditions);

				break;
			case MMC_QUAD:
				//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
				// Instead of:
				// mmr_fa_node_coor(meshId, faces[i], indices, coords);
				// quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 4);
				mmr_fa_node_coor(meshId, faces[i], indices, nullptr);

				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				mmr_node_coor(meshId, indices[2], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				mmr_node_coor(meshId, indices[3], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				mmr_node_coor(meshId, indices[4], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadCoords);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);

				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->lineCoords);
				//</ModFEM.Mesh-1.workaround>

				appendRealVector<INDEX_DIM * 4>(indices + 1, m->quadIndices);	// First element contains amount of coordinates.

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, vector3);
				for (int i = 0; i < 4; i++)
					appendRealVector<NORMAL_DIM>(vector3, m->quadNormals);

				bc = mmr_fa_bc(meshId, faces[i]);
				for (int i = 0; i < 4; i++)
					appendScalar<BOUNDARY_CONDITION_DIM>(bc, m->quadBoundaryConditions);


				// Tessellation of quads due to lack of GL_QUADS in modern OpenGL.

				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 1, m->quadTriangleIndices);

				mmr_node_coor(meshId, indices[2], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 2, m->quadTriangleIndices);

				mmr_node_coor(meshId, indices[3], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 3, m->quadTriangleIndices);

				mmr_node_coor(meshId, indices[4], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 4, m->quadTriangleIndices);

				mmr_node_coor(meshId, indices[3], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 3, m->quadTriangleIndices);

				mmr_node_coor(meshId, indices[1], vector3);
				appendRealVector<COORD_DIM>(vector3, m->quadTriangleCoords);
				appendRealVector<INDEX_DIM>(indices + 1, m->quadTriangleIndices);

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, vector3);
				for (int i = 0; i < 6; i++)
					appendRealVector<NORMAL_DIM>(vector3, m->quadTriangleNormals);

				bc = mmr_fa_bc(meshId, faces[i]);
				for (int i = 0; i < 6; i++)
					appendScalar<BOUNDARY_CONDITION_DIM>(bc, m->quadTriangleBoundaryConditions);

				break;
			default:
				CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faces[i]) << "'");
		}
	}
}

void ElementData::updateTriangleArrays(int meshId, int faceIndex)
{
	double vector3[3];	// Space required for x, y, z coordinates.
	double area;		// Face area.
	int indices[4];		// Triangle faces require 4 indices (first element contains amount of coordinates).
	int bc;				// Boundary condition flag.

	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// mmr_fa_node_coor(meshId, faceId, indices, coords);
	// triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 3);
	mmr_fa_node_coor(meshId, faceIndex, indices, nullptr);

	mmr_node_coor(meshId, indices[1], vector3);
	appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendRealVector<COORD_DIM>(vector3, m->lineCoords);

	mmr_node_coor(meshId, indices[2], vector3);
	appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[3], vector3);
	appendRealVector<COORD_DIM>(vector3, m->triangleCoords);
	appendRealVector<COORD_DIM>(vector3, m->lineCoords);

	appendRealVector<COORD_DIM>(vector3, m->lineCoords);
	mmr_node_coor(meshId, indices[1], vector3);
	appendRealVector<COORD_DIM>(vector3, m->lineCoords);
	//</ModFEM.Mesh-1.workaround>

	appendRealVector<INDEX_DIM * 3>(indices + 1, m->triangleIndices);	// First element contains amount of coordinates.

	// Append normal to each vertex.
	mmr_fa_area(meshId, faceIndex, & area, vector3);
	for (int i = 0; i < 3; i++)
		appendRealVector<NORMAL_DIM>(vector3, m->triangleNormals);

	// Append boundary condition flags to each vertex.

	bc = mmr_fa_bc(meshId, faceIndex);
	for (int i = 0; i < 3; i++)
		appendScalar<BOUNDARY_CONDITION_DIM>(bc, m->triangleBoundaryConditions);
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

}
}
