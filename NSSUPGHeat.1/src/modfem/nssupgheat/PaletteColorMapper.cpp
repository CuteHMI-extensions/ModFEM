#include "../../../include/modfem/nssupgheat/PaletteColorMapper.hpp"

namespace modfem {
namespace nssupgheat {

const QColor PaletteColorMapper::INITIAL_DEFAULT_COLOR = QColor("grey");

PaletteColorMapper::PaletteColorMapper(QObject * parent):
	AbstractColorMapper(parent),
	m(new Members{
	INITIAL_DEFAULT_COLOR,
	{}})
{
	connect(this, & AbstractColorMapper::inputTypeChanged, this, & PaletteColorMapper::updateOutput);
	connect(this, & AbstractColorMapper::inputChanged, this, & PaletteColorMapper::updateOutput);
	connect(this, & PaletteColorMapper::paletteChanged, this, & PaletteColorMapper::updateOutput);
}

QColor PaletteColorMapper::defaultColor() const
{
	return m->defaultColor;
}

void PaletteColorMapper::setDefaultColor(QColor defaultColor)
{
	if (m->defaultColor != defaultColor) {
		m->defaultColor = defaultColor;
		emit defaultColorChanged();
	}
}

QVariantList PaletteColorMapper::palette() const
{
	return m->palette;
}

void PaletteColorMapper::setPalette(QVariantList map)
{
	if (m->palette != map) {
		m->palette = map;
		emit paletteChanged();
	}
}

void PaletteColorMapper::updateOutput()
{
	output().clear();
	output().reserve(count() * OUTPUT_SIZE);
	QByteArray::const_iterator indicesBytes = input().constBegin();
	while (indicesBytes != input().constEnd()) {
		int colorIndex = *reinterpret_cast<const int *>(indicesBytes);
		QColor color = m->palette.value(colorIndex, m->defaultColor).value<QColor>();
//		CUTEHMI_DEBUG("colorIndex: " << colorIndex << " " << color);
		float colorVector[3];
		colorVector[0] = color.redF();
		colorVector[1] = color.greenF();
		colorVector[2] = color.blueF();
		output().append(reinterpret_cast<char *>(colorVector), OUTPUT_SIZE);
		indicesBytes += inputTypeSize();
	}

	emit outputChanged();
}

}
}
