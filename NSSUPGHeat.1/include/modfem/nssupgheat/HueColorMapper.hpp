#ifndef HUECOLORMAPPER_HPP
#define HUECOLORMAPPER_HPP

#include "internal/common.hpp"

#include <QObject>
#include <QColor>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API HueColorMapper:
	public QObject
{
		Q_OBJECT

	public:
		enum InputType {
			INPUT_FLOAT,
			INPUT_DOUBLE
		};
		Q_ENUM(InputType)

		static constexpr InputType INITIAL_INPUT_TYPE = INPUT_FLOAT;

		static constexpr qreal INITIAL_SATURATION = 1.0;

		static constexpr qreal INITIAL_LIGHTNESS = 0.5;

		static constexpr qreal INITIAL_HUE_BEGIN = 240.0 / 360.0;

		static constexpr qreal INITIAL_HUE_END = 1.0;

		Q_PROPERTY(InputType inputType READ inputType WRITE setInputType NOTIFY inputTypeChanged)

		Q_PROPERTY(QByteArray input READ input WRITE setInput NOTIFY inputChanged)

		Q_PROPERTY(QByteArray output READ output NOTIFY outputChanged)

		Q_PROPERTY(qreal saturation READ saturation WRITE setSaturation NOTIFY saturationChanged)

		Q_PROPERTY(qreal lightness READ lightness WRITE setLightness NOTIFY lightnessChanged)

		Q_PROPERTY(qreal hueBegin READ hueBegin WRITE setHueBegin NOTIFY hueBeginChanged)

		Q_PROPERTY(qreal hueEnd READ hueEnd WRITE setHueEnd NOTIFY hueEndChanged)

		Q_PROPERTY(qreal valueBegin READ valueBegin WRITE setValueBegin NOTIFY valueBeginChanged)

		Q_PROPERTY(qreal valueEnd READ valueEnd WRITE setValueEnd NOTIFY valueEndChanged)

		Q_PROPERTY(int count READ count NOTIFY countChanged)

		HueColorMapper(QObject * parent = nullptr);

		InputType inputType() const;

		void setInputType(InputType type);

		QByteArray input() const;

		void setInput(const QByteArray & input);

		QByteArray output() const;

		qreal saturation() const;

		void setSaturation(qreal saturation);

		qreal lightness() const;

		void setLightness(qreal lightness);

		qreal hueBegin() const;

		void setHueBegin(qreal hueBegin);

		qreal hueEnd() const;

		void setHueEnd(qreal hueEnd);

		qreal valueBegin() const;

		void setValueBegin(qreal valueBegin);

		qreal valueEnd() const;

		void setValueEnd(qreal valueEnd);

		int count() const;

	signals:
		void inputTypeChanged();

		void inputChanged();

		void outputChanged();

		void saturationChanged();

		void lightnessChanged();

		void hueBeginChanged();

		void hueEndChanged();

		void valueBeginChanged();

		void valueEndChanged();

		void countChanged();

	private:
		static constexpr std::size_t COLOR_VECTOR_DIM = 3;

		static constexpr std::size_t COLOR_VECTOR_SIZE = sizeof(float) * COLOR_VECTOR_DIM;

		void setCount(int count);

		void updateOutput();

		std::size_t inputTypeSize() const;

		QColor pullColorFromInput(QByteArray::const_iterator & pos, QByteArray::const_iterator end);

		struct Members {
			InputType inputType;
			qreal saturation;
			qreal lightness;
			qreal hueBegin;
			qreal hueEnd;
			qreal valueBegin;
			qreal valueEnd;
			QByteArray input;
			QByteArray output;
			int count;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif
