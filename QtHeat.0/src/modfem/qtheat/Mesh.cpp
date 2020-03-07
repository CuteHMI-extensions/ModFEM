#include "../../../include/modfem/qtheat/Mesh.hpp"

#include <modfem/mmh_intf.h>

namespace modfem {
namespace qtheat {

Mesh::Mesh(QObject * parent):
	QObject(parent),
	m(new Members{
	0,
	{},
	new FaceData(this)})
{
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

QVariantMap Mesh::triangles() const
{
	return m->triangles;
}

void Mesh::init(int meshId)
{
	setNodeData(meshId);
	m->faceData->init(meshId);
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

}
}
