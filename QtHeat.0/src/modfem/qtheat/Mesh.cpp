#include "../../../include/modfem/qtheat/Mesh.hpp"

namespace modfem {
namespace qtheat {

Mesh::Mesh(QObject * parent):
	QObject(parent),
	m(new Members{
	0, {}})
{
}

int Mesh::nodeCount() const
{
	return m->nodeCount;
}

QByteArray Mesh::nodeData() const
{
	return m->nodeData;
}

void Mesh::init(int meshId)
{
	setNodeData(meshId);
}

void Mesh::setNodeData(int meshId)
{
	//temp
//	QDataStream stream(& m->meshData, QIODevice::WriteOnly);
//	stream.setFloatingPointPrecision(QDataStream::SinglePrecision);

//	stream << (qint32) - 5.0f << (qint32) - 5.0f << (qint32)0.0f;
//	stream << (qint32)0.0f << (qint32) - 5.0f << (qint32)0.0f;
//	stream << -5.0f << -5.0f << 0.0f;
//	m->buffer->setData(tempData);
	double point = -5.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = -5.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	point = 0.0;
	m->nodeData.append(reinterpret_cast<char *>(& point), sizeof (double));
	emit nodeDataChanged();
	setNodeCount(7);
	//endtemp
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
