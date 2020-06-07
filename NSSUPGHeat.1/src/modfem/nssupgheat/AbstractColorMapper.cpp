#include <modfem/nssupgheat/AbstractColorMapper.hpp>

namespace modfem {
namespace nssupgheat {

constexpr AbstractColorMapper::InputType AbstractColorMapper::INITIAL_INPUT_TYPE;
constexpr std::size_t AbstractColorMapper::OUTPUT_DIM;
constexpr std::size_t AbstractColorMapper::OUTPUT_SIZE;

AbstractColorMapper::AbstractColorMapper(QObject * parent):
	QObject(parent),
	m(new Members{
	INITIAL_INPUT_TYPE,
	{},
	{},
	0})
{
}

AbstractColorMapper::InputType AbstractColorMapper::inputType() const
{
	return m->inputType;
}

void AbstractColorMapper::setInputType(AbstractColorMapper::InputType type)
{
	if (m->inputType != type) {
		m->inputType = type;
		emit inputTypeChanged();
	}
}

QByteArray AbstractColorMapper::input() const
{
	return m->input;
}

void AbstractColorMapper::setInput(const QByteArray & input)
{
	if (m->input.constData() != input.constData()) {
		m->input = input;
		emit inputChanged();
		setCount(m->input.length() / inputTypeSize());
	}
}

QByteArray AbstractColorMapper::output() const
{
	return m->output;
}

int AbstractColorMapper::count() const
{
	return m->count;
}

void AbstractColorMapper::setCount(int count)
{
	if (m->count != count) {
		m->count = count;
		emit countChanged();
	}
}

std::size_t AbstractColorMapper::inputTypeSize() const
{
	switch (m->inputType) {
		case INPUT_FLOAT:
			return sizeof(float);
		case INPUT_DOUBLE:
			return sizeof(double);
		case INPUT_INT:
			return sizeof(int);
		default:
			CUTEHMI_CRITICAL("Unrecognized input type '" << m->inputType << "'.");
	}
	return sizeof(float);
}

QByteArray & AbstractColorMapper::input()
{
	return m->input;
}

QByteArray & AbstractColorMapper::output()
{
	return m->output;
}

}
}
