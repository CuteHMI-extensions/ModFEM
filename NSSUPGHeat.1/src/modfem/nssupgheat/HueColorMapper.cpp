#include "../../../include/modfem/nssupgheat/HueColorMapper.hpp"

namespace modfem {
namespace nssupgheat {

constexpr qreal HueColorMapper::INITIAL_HUE_BEGIN;
constexpr qreal HueColorMapper::INITIAL_HUE_END;
constexpr qreal HueColorMapper::INITIAL_SATURATION;
constexpr qreal HueColorMapper::INITIAL_LIGHTNESS;

HueColorMapper::HueColorMapper(QObject * parent):
	AbstractColorMapper(parent),
	m(new Members{
	INITIAL_SATURATION,
	INITIAL_LIGHTNESS,
	INITIAL_HUE_BEGIN,
	INITIAL_HUE_END,
	0.0,
	0.0})
{
	connect(this, & AbstractColorMapper::inputTypeChanged, this, & HueColorMapper::updateOutput);
	connect(this, & AbstractColorMapper::inputChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::saturationChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::lightnessChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::hueBeginChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::hueEndChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::valueBeginChanged, this, & HueColorMapper::updateOutput);
	connect(this, & HueColorMapper::valueEndChanged, this, & HueColorMapper::updateOutput);
}

qreal HueColorMapper::saturation() const
{
	return m->saturation;
}

void HueColorMapper::setSaturation(qreal saturation)
{
	if (m->saturation != saturation) {
		m->saturation = saturation;
		emit saturationChanged();
	}
}

qreal HueColorMapper::lightness() const
{
	return m->lightness;
}

void HueColorMapper::setLightness(qreal lightness)
{
	if (m->lightness != lightness) {
		m->lightness = lightness;
		emit lightnessChanged();
	}
}

qreal HueColorMapper::hueBegin() const
{
	return m->hueBegin;
}

void HueColorMapper::setHueBegin(qreal fromHue)
{
	if (m->hueBegin != fromHue) {
		m->hueBegin = fromHue;
		emit hueBeginChanged();
	}
}

qreal HueColorMapper::hueEnd() const
{
	return m->hueEnd;
}

void HueColorMapper::setHueEnd(qreal toHue)
{
	if (m->hueEnd != toHue) {
		m->hueEnd = toHue;
		emit hueEndChanged();
	}
}

qreal HueColorMapper::valueBegin() const
{
	return m->valueBegin;
}

void HueColorMapper::setValueBegin(qreal from)
{
	if (m->valueBegin != from) {
		m->valueBegin = from;
		emit valueBeginChanged();
	}
}

qreal HueColorMapper::valueEnd() const
{
	return m->valueEnd;
}

void HueColorMapper::setValueEnd(qreal to)
{
	if (m->valueEnd != to) {
		m->valueEnd = to;
		emit valueEndChanged();
	}
}

void HueColorMapper::updateOutput()
{
	output().clear();
	output().reserve(count() * OUTPUT_SIZE);
	QColor color;
	QByteArray::const_iterator inputIterator = input().constBegin();
	while ((color = pullColorFromInput(inputIterator, input().constEnd())).isValid()) {
		float colorVector[3];
//		color = QColor::fromRgb(0, 120, 233);	//TEMP
//		CUTEHMI_DEBUG("color; " << color);
		colorVector[0] = color.redF();
		colorVector[1] = color.greenF();
		colorVector[2] = color.blueF();
		output().append(reinterpret_cast<char *>(colorVector), OUTPUT_SIZE);
	}
	emit outputChanged();
}

QColor HueColorMapper::pullColorFromInput(QByteArray::const_iterator & pos, QByteArray::const_iterator end)
{
	qreal value;

	if (pos != end) {
		switch (inputType()) {
			case INPUT_FLOAT: {
				value = *reinterpret_cast<const float *>(pos);
				pos += sizeof(float);
				break;
			}
			case INPUT_DOUBLE: {
				value = *reinterpret_cast<const double *>(pos);
				pos += sizeof(double);
				break;
			}
			case INPUT_INT: {
				value = *reinterpret_cast<const int *>(pos);
				pos += sizeof(int);
				break;
			}
			default:
				CUTEHMI_CRITICAL("Unrecognized input type '" << inputType() << "'.");
				return QColor();
		}
		qreal valuePercent = (value - m->valueBegin) / (m->valueEnd - m->valueBegin);
		qreal hue = m->hueBegin + valuePercent * (m->hueEnd - m->hueBegin);
//		CUTEHMI_DEBUG("Mapping value " << value << " value percentage: " << valuePercent << " to hue: " << hue);
		return QColor::fromHslF(qMax(hueBegin(), qMin(hue, hueEnd())), saturation(), lightness());
	}
	return QColor();
}

}
}
