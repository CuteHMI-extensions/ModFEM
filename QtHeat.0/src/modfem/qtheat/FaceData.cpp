#include "../../../include/modfem/qtheat/FaceData.hpp"

#include <modfem/mmh_intf.h>

namespace modfem {
namespace qtheat {

FaceData::FaceData(QObject * parent):
	QObject(parent),
	m(new Members{
	0,
	{},
	{},
	{},
	0,
	{},
	{}})
{
}

int FaceData::triangleCount() const
{
	return m->triangleCount;
}

QByteArray FaceData::triangleIndices() const
{
	return m->triangleIndices;
}

QByteArray FaceData::triangleCoords() const
{
	return m->triangleCoords;
}

QByteArray FaceData::triangleNormals() const
{
	return m->triangleNormals;
}

int FaceData::quadCount() const
{
	return m->quadCount;
}

QByteArray FaceData::quadIndices() const
{
	return m->quadIndices;
}

QByteArray FaceData::quadCoords()
{
	return m->quadCoords;
}

void FaceData::init(int meshId)
{
	m->triangleCoords.clear();
	m->triangleIndices.clear();
	m->triangleNormals.clear();
	m->quadCoords.clear();
	m->quadIndices.clear();

	// First pass - calculate number of faces of various types.
	int triangleCount = 0;
	int quadCount = 0;
	int otherCount = 0;
	int faceId = mmr_get_next_act_face(meshId, 0);
	while (faceId != 0) {
		switch (mmr_fa_type(meshId, faceId)) {
			case MMC_TRIA:
				triangleCount++;
				break;
			case MMC_QUAD:
				quadCount++;
				break;
			default:
				otherCount++;
				CUTEHMI_CRITICAL("Unsupported face type: '" << mmr_fa_type(meshId, faceId) << "'");
		}
		faceId = mmr_get_next_act_face(meshId, faceId);
	}
	CUTEHMI_ASSERT(mmr_get_nr_face(meshId) == triangleCount + quadCount + otherCount,
			QString("Amount of active faces mismatch (1st method gives %1, 2nd method gives %2)").arg(mmr_get_nr_face(meshId)).arg(mmr_get_nr_face(meshId)).toLocal8Bit().constData());

	// Second pass - fill in arrays.
	m->triangleCoords.reserve(sizeof(double) * triangleCount * 3);
	m->triangleIndices.reserve(sizeof(int) * triangleCount * 3);
	m->triangleNormals.reserve(sizeof(double) * triangleCount * 3);
	m->quadCoords.reserve(sizeof(double) * quadCount * 4);
	m->quadIndices.reserve(sizeof(int) * quadCount * 4);
	double coords[4];	// Quad faces require 4 coordinates.
	int indices[5];		// Quad faces require 5 indices (first element contains amount of coordinates).
	double normal[3];	// Normal vector, perpendicular to surface.
	double area;		// Face area.
	faceId = mmr_get_next_act_face(meshId, 0);
	while (faceId != 0) {
		switch (mmr_fa_type(meshId, faceId)) {
			case MMC_TRIA:
				//<workaround target="ModFEM-mm_t4_prism" cause="bug">
				// Instead of:
				// mmr_fa_node_coor(meshId, faceId, indices, coords);
				mmr_fa_node_coor(meshId, faceId, indices, nullptr);

//				CUTEHMI_DEBUG("Appending " << faceId << " coords: " << coords[0] << " " << coords[1] << " " << coords[2]);
				// workaround
				for (int i = 1; i < 4; i++) {
					mmr_node_coor(meshId, indices[i], coords);
//					CUTEHMI_DEBUG("Appending " << faceId << " coords: " << coords[0] << " " << coords[1] << " " << coords[2]);
//					CUTEHMI_DEBUG("Appending " << faceId << " coords: " << coords[0] << " " << coords[1] << " " << coords[2]);
					m->triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
				}

//				m->triangleCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 3);
//				CUTEHMI_DEBUG("Appending " << faceId << " indices: " << indices[0] << " " << indices[1] << " " << indices[2] << " " << indices[3]);
				m->triangleIndices.append(reinterpret_cast<char *>(indices + 1), sizeof(int) * 3);

				mmr_fa_area(meshId, faceId, & area, normal);
				// Append normal for each vertex.
				for (int i = 0; i < 3; i++)
					m->triangleNormals.append(reinterpret_cast<char *>(normal), sizeof(double) * 3);

				break;
			case MMC_QUAD:
				mmr_fa_node_coor(meshId, faceId, indices, coords);
				m->quadCoords.append(reinterpret_cast<char *>(coords), sizeof(double) * 4);
				m->quadIndices.append(reinterpret_cast<char *>(indices), sizeof(int) * 4);
				break;
			default:
				CUTEHMI_WARNING("Unsupported face type: '" << mmr_fa_type(meshId, faceId) << "'");
		}
		faceId = mmr_get_next_act_face(meshId, faceId);
	}

	setTriangleCount(triangleCount);
	emit triangleIndicesChanged();
	emit triangleCoordsChanged();
	emit triangleNormalsChanged();
	setQuadCount(quadCount);
	emit quadIndicesChanged();
	emit quadCoordsChanged();
}

void FaceData::setTriangleCount(int count)
{
	if (m->triangleCount != count) {
		m->triangleCount = count;
		emit triangleCountChanged();
	}
}

void FaceData::setQuadCount(int count)
{
	if (m->quadCount != count) {
		m->quadCount = count;
		emit quadCountChanged();
	}
}

}
}
