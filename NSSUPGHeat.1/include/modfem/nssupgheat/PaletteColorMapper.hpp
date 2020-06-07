#ifndef VERTEXCOLORMAPPER_HPP
#define VERTEXCOLORMAPPER_HPP

#include "AbstractColorMapper.hpp"

#include <QColor>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API PaletteColorMapper:
	public AbstractColorMapper
{
		Q_OBJECT

	public:
		static const QColor INITIAL_DEFAULT_COLOR;

		Q_PROPERTY(QColor defaultColor READ defaultColor WRITE setDefaultColor NOTIFY defaultColorChanged)

		Q_PROPERTY(QVariantList palette READ palette WRITE setPalette NOTIFY paletteChanged)

		PaletteColorMapper(QObject * parent = nullptr);

		QColor defaultColor() const;

		void setDefaultColor(QColor defaultColor);

		QVariantList palette() const;

		void setPalette(QVariantList palette);

	signals:
		void defaultColorChanged();

		void paletteChanged();

	private:
		void updateOutput();

		struct Members {
			QColor defaultColor;
			QVariantList palette;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // VERTEXCOLORMAPPER_HPP
