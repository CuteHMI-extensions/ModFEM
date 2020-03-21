#include "../../../include/modfem/heat/ElementData.hpp"

#include <modfem/mmh_intf.h>

namespace modfem {
namespace heat {

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
	clearArrays();
	reserveArrays();
	updateArrays(meshId);
}

void ElementData::clearArrays()
{
	m->lineCoords.clear();
	m->triangleCoords.clear();
	m->triangleIndices.clear();
	m->triangleNormals.clear();
	m->quadCoords.clear();
	m->quadIndices.clear();
	m->quadNormals.clear();
	m->quadTriangleCoords.clear();
	m->quadTriangleIndices.clear();
	m->quadTriangleNormals.clear();
}

void ElementData::reserveArrays()
{
	m->lineCoords.reserve(m->count["lines"].toInt() * 2 * COORD_SIZE);
	m->triangleCoords.reserve(m->count["triangles"].toInt() * 3 * COORD_SIZE);
	m->triangleIndices.reserve(m->count["triangles"].toInt() * 3 * INDEX_SIZE);
	m->triangleNormals.reserve(m->count["triangles"].toInt() * 3 * NORMAL_SIZE);
	m->quadCoords.reserve(m->count["quads"].toInt() * 4 * COORD_SIZE);
	m->quadIndices.reserve(m->count["quads"].toInt() * 4 * INDEX_SIZE);
	m->quadNormals.reserve(m->count["quads"].toInt() * 4 * NORMAL_SIZE);
	m->quadTriangleCoords.reserve(m->count["quads"].toInt() * 6 * COORD_SIZE);
	m->quadTriangleIndices.reserve(m->count["quads"].toInt() * 6 * INDEX_SIZE);
	m->quadTriangleNormals.reserve(m->count["quads"].toInt() * 6 * NORMAL_SIZE);
	m->nodeCoords.reserve(m->nodes["count"].toInt() * COORD_SIZE);	/// @todo adding a nodes per face requires different value.
}

void ElementData::updateArrays(int meshId)
{
	int elementId = mmr_get_next_act_elem(meshId, 0);

	while (elementId != 0) {
		updateArrays(meshId, elementId);
		elementId = mmr_get_next_act_elem(meshId, elementId);
	}

	/// @todo add nodes per face.
	double coords[3];
	for (int i = 1; i <= mmr_get_max_node_id(meshId); i++) {
		if (mmr_node_status(meshId, i) == MMC_ACTIVE) {
			mmr_node_coor(meshId, i, coords);
			m->nodeCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
		}
	}

	updateProperties();
}

void ElementData::updateArrays(int meshId, int elementId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces in 1st element).

	int faces[FACES_MAX] = {0};
	int orientation[FACES_MAX];
	//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
	// Instead of:
	// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
	double coords[3];	// Space required for x, y, z coordinates.
	//</ModFEM.Mesh-1.workaround>
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	double normal[3];	// Normal vector, perpendicular to surface.
	double area;		// Face area.

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

				mmr_node_coor(meshId, indices[1], coords);
				m->triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				mmr_node_coor(meshId, indices[2], coords);
				m->triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				mmr_node_coor(meshId, indices[3], coords);
				m->triangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				mmr_node_coor(meshId, indices[1], coords);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				//</ModFEM.Mesh-1.workaround>

				m->triangleIndices.append(reinterpret_cast<char *>(indices + 1), INDEX_SIZE * 3);	// First element contains amount of coordinates.

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, normal);
				for (int i = 0; i < 3; i++)
					m->triangleNormals.append(reinterpret_cast<char *>(normal), NORMAL_SIZE);

				break;
			case MMC_QUAD:
				//<ModFEM.Mesh-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
				// Instead of:
				// mmr_fa_node_coor(meshId, faces[i], indices, coords);
				// quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE * 4);
				mmr_fa_node_coor(meshId, faces[i], indices, nullptr);

				mmr_node_coor(meshId, indices[1], coords);
				m->quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				mmr_node_coor(meshId, indices[2], coords);
				m->quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				mmr_node_coor(meshId, indices[3], coords);
				m->quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				mmr_node_coor(meshId, indices[4], coords);
				m->quadCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);

				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				mmr_node_coor(meshId, indices[1], coords);
				m->lineCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				//</ModFEM.Mesh-1.workaround>

				m->quadIndices.append(reinterpret_cast<char *>(indices + 1), INDEX_SIZE * 4);	// First element contains amount of coordinates.

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, normal);
				for (int i = 0; i < 4; i++)
					m->quadNormals.append(reinterpret_cast<char *>(normal), NORMAL_SIZE);


				// Tessellation of quads due to lack of GL_QUADS in modern OpenGL.

				mmr_node_coor(meshId, indices[1], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 1), INDEX_SIZE);

				mmr_node_coor(meshId, indices[2], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 2), INDEX_SIZE);

				mmr_node_coor(meshId, indices[3], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 3), INDEX_SIZE);

				mmr_node_coor(meshId, indices[4], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 4), INDEX_SIZE);

				mmr_node_coor(meshId, indices[3], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 3), INDEX_SIZE);

				mmr_node_coor(meshId, indices[1], coords);
				m->quadTriangleCoords.append(reinterpret_cast<char *>(coords), COORD_SIZE);
				m->quadTriangleIndices.append(reinterpret_cast<char *>(indices + 1), INDEX_SIZE);

				// Append normal to each vertex.
				mmr_fa_area(meshId, faces[i], & area, normal);
				for (int i = 0; i < 6; i++)
					m->quadTriangleNormals.append(reinterpret_cast<char *>(normal), NORMAL_SIZE);

				break;
			default:
				CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faces[i]) << "'");
		}
	}
}

void ElementData::updateProperties()
{
	m->triangles["coords"] = m->triangleCoords;
	m->triangles["count"] = m->triangleCoords.count() / (COORD_SIZE * 3);
	m->triangles["normals"] = m->triangleNormals;
	m->quads["coords"] = m->quadCoords;
	m->quads["count"] = m->quadCoords.count() / (COORD_SIZE * 4);
	m->quads["normals"] = m->quadNormals;
	m->quads["triangleCoords"] = m->quadTriangleCoords;
	m->quads["triangleCount"] = m->quadTriangleCoords.count() / (COORD_SIZE * 3);
	m->quads["triangleNormals"] = m->quadTriangleNormals;
	m->lines["coords"] = m->lineCoords;
	m->lines["count"] = m->lineCoords.count() / (COORD_SIZE * 2);
	m->nodes["coords"] = m->nodeCoords;
	m->nodes["count"] = m->nodeCoords.count() / (COORD_SIZE);

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
