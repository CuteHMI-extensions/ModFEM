#include "../../../include/modfem/nssupgheat/HueColorMapper.hpp"

namespace modfem {
namespace nssupgheat {

constexpr HueColorMapper::InputType HueColorMapper::INITIAL_INPUT_TYPE;
constexpr qreal HueColorMapper::INITIAL_HUE_BEGIN;
constexpr qreal HueColorMapper::INITIAL_HUE_END;
constexpr qreal HueColorMapper::INITIAL_SATURATION;
constexpr qreal HueColorMapper::INITIAL_LIGHTNESS;

HueColorMapper::HueColorMapper(QObject * parent):
	QObject(parent),
	m(new Members{
	INITIAL_INPUT_TYPE,
	INITIAL_SATURATION,
	INITIAL_LIGHTNESS,
	INITIAL_HUE_BEGIN,
	INITIAL_HUE_END,
	0.0,
	0.0,
	{},
	{},
	0})
{
}

HueColorMapper::InputType HueColorMapper::inputType() const
{
	return m->inputType;
}

void HueColorMapper::setInputType(HueColorMapper::InputType type)
{
	if (m->inputType != type) {
		m->inputType = type;
		emit inputTypeChanged();
		updateOutput();
	}
}

QByteArray HueColorMapper::input() const
{
	return m->input;
}

void HueColorMapper::setInput(const QByteArray & input)
{
	if (m->input.constData() != input.constData()) {
		m->input = input;
		emit inputChanged();
		setCount(m->input.length() / inputTypeSize());
		updateOutput();
	}
}

QByteArray HueColorMapper::output() const
{
	return m->output;
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
		updateOutput();
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
		updateOutput();
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
		updateOutput();
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
		updateOutput();
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
		updateOutput();
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
		updateOutput();
	}
}

int HueColorMapper::count() const
{
	return m->count;
}

void HueColorMapper::setCount(int count)
{
	if (m->count != count) {
		m->count = count;
		emit countChanged();
	}
}

void HueColorMapper::updateOutput()
{
	m->output.clear();
	m->output.reserve(count() * COLOR_VECTOR_SIZE);
	QColor color;
	QByteArray::const_iterator inputIterator = m->input.constBegin();
	while ((color = pullColorFromInput(inputIterator, m->input.constEnd())).isValid()) {
		float colorVector[3];
//		color = QColor::fromRgb(0, 120, 233);	//TEMP
//		CUTEHMI_DEBUG("color; " << color);
		colorVector[0] = color.redF();
		colorVector[1] = color.greenF();
		colorVector[2] = color.blueF();
		m->output.append(reinterpret_cast<char *>(colorVector), COLOR_VECTOR_SIZE);
	}
	emit outputChanged();
}

std::size_t HueColorMapper::inputTypeSize() const
{
	switch (m->inputType) {
		case INPUT_FLOAT:
			return sizeof(float);
		case INPUT_DOUBLE:
			return sizeof(double);
		default:
			CUTEHMI_CRITICAL("Unrecognized input type '" << m->inputType << "'.");
	}
	return sizeof(float);
}

QColor HueColorMapper::pullColorFromInput(QByteArray::const_iterator & pos, QByteArray::const_iterator end)
{
	qreal value;

	if (pos != end) {
		switch (m->inputType) {
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
			default:
				CUTEHMI_CRITICAL("Unrecognized input type '" << m->inputType << "'.");
				return QColor();
		}
		qreal valuePercent = (value - m->valueBegin) / (m->valueEnd - m->valueBegin);
		qreal hue = m->hueBegin + valuePercent * (m->hueEnd - m->hueBegin);
		CUTEHMI_DEBUG("Mapping value " << value << " value percentage: " << valuePercent << " to hue: " << hue);
		return QColor::fromHslF(qMax(hueBegin(), qMin(hue, hueEnd())), saturation(), lightness());
	}
	return QColor();
}

}
}
