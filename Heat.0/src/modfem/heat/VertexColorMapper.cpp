#include "../../../include/modfem/heat/VertexColorMapper.hpp"

namespace modfem {
namespace heat {

const QColor VertexColorMapper::INITIAL_DEFAULT_COLOR = QColor("limegreen");

VertexColorMapper::VertexColorMapper(QObject * parent):
	QObject(parent),
	m(new Members{
	INITIAL_DEFAULT_COLOR,
	{},
	{},
	{},
	0})
{
}

QColor VertexColorMapper::defaultColor() const
{
	return m->defaultColor;
}

void VertexColorMapper::setDefaultColor(QColor defaultColor)
{
	if (m->defaultColor != defaultColor) {
		m->defaultColor = defaultColor;
		emit defaultColorChanged();
	}
}

QVariantList VertexColorMapper::map() const
{
	return m->map;
}

void VertexColorMapper::setMap(QVariantList map)
{
	if (m->map != map) {
		m->map = map;
		emit mapChanged();
		updateColorVertices();
	}
}

QByteArray VertexColorMapper::colorIndices() const
{
	return m->colorIndices;
}

void VertexColorMapper::setColorIndices(const QByteArray & colorIndices)
{
	if (m->colorIndices.constData() != colorIndices.constData()) {
		m->colorIndices = colorIndices;
		emit colorIndicesChanged();
		setCount(m->colorIndices.length() / COLOR_INDEX_SIZE);
		updateColorVertices();
	}
}

QByteArray VertexColorMapper::colorVertices() const
{
	return m->colorVertices;
}

int VertexColorMapper::count() const
{
	return m->count;
}

void VertexColorMapper::setCount(int count)
{
	if (m->count != count) {
		m->count = count;
		emit countChanged();
	}
}

void VertexColorMapper::updateColorVertices()
{
	m->colorVertices.clear();
	m->colorVertices.reserve(count() * COLOR_VECTOR_SIZE);
	QByteArray::const_iterator indicesBytes = m->colorIndices.constBegin();
	while (indicesBytes != m->colorIndices.constEnd()) {
		int colorIndex = *reinterpret_cast<const int *>(indicesBytes);
		QColor color = m->map.value(colorIndex, m->defaultColor).value<QColor>();
		float colorVector[3];
		colorVector[0] = color.redF();
		colorVector[1] = color.greenF();
		colorVector[2] = color.blueF();
		m->colorVertices.append(reinterpret_cast<char *>(colorVector), COLOR_VECTOR_SIZE);
		indicesBytes += COLOR_INDEX_SIZE;
	}

	emit colorVerticesChanged();
}

}
}
