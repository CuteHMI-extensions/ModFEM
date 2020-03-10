#include "../../../include/modfem/qtheat/Mesh.hpp"

#include <modfem/mmh_intf.h>

namespace modfem {
namespace qtheat {

Mesh::Mesh(QObject * parent):
	QObject(parent),
	m(new Members{
	0,
	{},
	new FaceData(this),
	{},
	{},
	{},
	{}})
{
	m->triangles["count"] = 0;
	m->quads["count"] = 0;
	m->lines["count"] = 0;
	m->nodes["count"] = 0;
}

int Mesh::nodeCount() const
{
	return m->nodeCount;
}

QByteArray Mesh::nodeCoords() const
{
	return m->nodeCoords;
}

FaceData * Mesh::faceData() const
{
	return m->faceData;
}

QVariantMap Mesh::nodes() const
{
	return m->nodes;
}

QVariantMap Mesh::triangles() const
{
	return m->triangles;
}

QVariantMap Mesh::quads() const
{
	return m->quads;
}

QVariantMap Mesh::lines() const
{
	return m->lines;
}

void Mesh::init(int meshId)
{
	setNodeData(meshId);
	m->faceData->init(meshId);
	update(meshId);
}

void Mesh::setNodeData(int meshId)
{
	int activeNodes = mmr_get_nr_node(meshId);
	m->nodeCoords.clear();
	m->nodeCoords.reserve(sizeof(double) * activeNodes * 3);

	double coords[3];
	for (int i = 0; i <= mmr_get_max_node_id(meshId); i++) {
		if (mmr_node_status(meshId, i) == MMC_ACTIVE) {
			mmr_node_coor(meshId, i, coords);
			m->nodeCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
		}
	}
	setNodeCount(activeNodes);
	emit nodeCoordsChanged();
}

void Mesh::setNodeCount(int count)
{
	if (m->nodeCount != count) {
		m->nodeCount = count;
		emit nodeCountChanged();
	}
}

void Mesh::update(int meshId)
{
	static constexpr int FACES_MAX = 7;	// Assuming brick is largest element in terms of amount of faces and it has 6 faces (+1, because mmr_el_faces() stores number of faces in 1st element).

	count(meshId);

	QByteArray lineCoords;
	QByteArray triangleCoords;
	QByteArray triangleIndices;
	QByteArray triangleNormals;
	QByteArray quadCoords;
	QByteArray quadIndices;
	QByteArray quadNormals;

	lineCoords.reserve(sizeof(double) * (m->triangles["count"].toInt() * 6 + m->quads["count"].toInt() * 8));
	triangleCoords.reserve(sizeof(double) * m->triangles["count"].toInt() * 3);
	triangleIndices.reserve(sizeof(double) * m->triangles["count"].toInt() * 3);
	triangleNormals.reserve(sizeof(double) * m->triangles["count"].toInt() * 3);
	quadCoords.reserve(sizeof(double) * m->quads["count"].toInt() * 4);
	quadIndices.reserve(sizeof(double) * m->quads["count"].toInt() * 4);
	quadNormals.reserve(sizeof(double) * m->quads["count"].toInt() * 4);

	int elementId = mmr_get_next_act_elem(meshId, 0);

	while (elementId != 0) {

		int faces[FACES_MAX] = {0};
		int orientation[FACES_MAX];
		//<ModFEM.QtHeat-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
		// Instead of:
		// double coords[4 * 3];	// Quad faces require 4 coordinates (times x, y, z).
		double coords[3];	// Space required for x, y, z coordinates.
		//</ModFEM.QtHeat-1.workaround>
		int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
		double normal[3];	// Normal vector, perpendicular to surface.
		double area;		// Face area.

		mmr_el_faces(meshId, elementId, faces, orientation);
		CUTEHMI_ASSERT(faces[0] < FACES_MAX, "not enough space reserved for 'faces' array");

		int faceCount = faces[0];
		for (int i = 1; i <= faceCount; i++) {
			switch (mmr_fa_type(meshId, faces[i])) {
				case MMC_TRIA:
					//<ModFEM.QtHeat-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
					// Instead of:
					// mmr_fa_node_coor(meshId, faceId, indices, coords);
					// triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3 * 3);
					mmr_fa_node_coor(meshId, faces[i], indices, nullptr);

					mmr_node_coor(meshId, indices[1], coords);
					triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					mmr_node_coor(meshId, indices[2], coords);
					triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					mmr_node_coor(meshId, indices[3], coords);
					triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					mmr_node_coor(meshId, indices[1], coords);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					//</ModFEM.QtHeat-1.workaround>

					triangleIndices.append(reinterpret_cast<char *>(indices + 1), sizeof(int) * 3);	// First element contains amount of coordinates.

					// Append normal to each vertex.
					mmr_fa_area(meshId, faces[i], & area, normal);
					for (int i = 0; i < 3; i++)
						triangleNormals.append(reinterpret_cast<char *>(normal), sizeof(double) * 3);

					break;
				case MMC_QUAD:
					//<ModFEM.QtHeat-1.workaround target="ModFEM-mm_t4_prism" cause="bug">
					// Instead of:
					// mmr_fa_node_coor(meshId, faceId, indices, coords);
					// quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 4 * 3);
					mmr_fa_node_coor(meshId, faces[i], indices, nullptr);

					mmr_node_coor(meshId, indices[1], coords);
					quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					mmr_node_coor(meshId, indices[2], coords);
					quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					mmr_node_coor(meshId, indices[3], coords);
					quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					mmr_node_coor(meshId, indices[4], coords);
					quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);

					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					mmr_node_coor(meshId, indices[1], coords);
					lineCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
					//</ModFEM.QtHeat-1.workaround>

					quadIndices.append(reinterpret_cast<char *>(indices + 1), sizeof(int) * 4);	// First element contains amount of coordinates.

					// Append normal to each vertex.
					mmr_fa_area(meshId, faces[i], & area, normal);
					for (int i = 0; i < 4; i++)
						quadNormals.append(reinterpret_cast<char *>(normal), sizeof(double) * 3);

					break;
				default:
					CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faces[i]) << "'");
			}
		}

		elementId = mmr_get_next_act_elem(meshId, elementId);
	}

	QByteArray nodeCoords;
	nodeCoords.reserve(sizeof(double) * m->nodes["count"].toInt() * 3);
	double coords[3];
	for (int i = 0; i <= mmr_get_max_node_id(meshId); i++) {
		if (mmr_node_status(meshId, i) == MMC_ACTIVE) {
			mmr_node_coor(meshId, i, coords);
			nodeCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
		}
	}

	m->triangles["coords"] = triangleCoords;
	m->triangles["normals"] = triangleNormals;
	m->quads["coords"] = quadCoords;
	m->quads["normals"] = quadNormals;
	m->lines["coords"] = lineCoords;
	m->nodes["coords"] = nodeCoords;

	emit trianglesChanged();
	emit quadsChanged();
	emit linesChanged();
	emit nodesChanged();
}

void Mesh::count(int meshId)
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

	m->triangles["count"] = triangleCount;
	m->quads["count"] = quadCount;
	m->lines["count"] = lineCount;
	m->nodes["count"] = mmr_get_nr_node(meshId);
}

}
}
