#ifndef HUECOLORMAPPER_HPP
#define HUECOLORMAPPER_HPP

#include "AbstractColorMapper.hpp"

#include <QObject>
#include <QColor>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API HueColorMapper:
	public AbstractColorMapper
{
		Q_OBJECT

	public:
		static constexpr qreal INITIAL_SATURATION = 1.0;

		static constexpr qreal INITIAL_LIGHTNESS = 0.5;

		static constexpr qreal INITIAL_HUE_BEGIN = 240.0 / 360.0;

		static constexpr qreal INITIAL_HUE_END = 1.0;

		Q_PROPERTY(qreal saturation READ saturation WRITE setSaturation NOTIFY saturationChanged)

		Q_PROPERTY(qreal lightness READ lightness WRITE setLightness NOTIFY lightnessChanged)

		Q_PROPERTY(qreal hueBegin READ hueBegin WRITE setHueBegin NOTIFY hueBeginChanged)

		Q_PROPERTY(qreal hueEnd READ hueEnd WRITE setHueEnd NOTIFY hueEndChanged)

		Q_PROPERTY(qreal valueBegin READ valueBegin WRITE setValueBegin NOTIFY valueBeginChanged)

		Q_PROPERTY(qreal valueEnd READ valueEnd WRITE setValueEnd NOTIFY valueEndChanged)

		HueColorMapper(QObject * parent = nullptr);

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

	signals:
		void saturationChanged();

		void lightnessChanged();

		void hueBeginChanged();

		void hueEndChanged();

		void valueBeginChanged();

		void valueEndChanged();

	private:
		void updateOutput();

		QColor pullColorFromInput(QByteArray::const_iterator & pos, QByteArray::const_iterator end);

		struct Members {
			qreal saturation;
			qreal lightness;
			qreal hueBegin;
			qreal hueEnd;
			qreal valueBegin;
			qreal valueEnd;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif
