#ifndef VERTEXCOLORMAPPER_HPP
#define VERTEXCOLORMAPPER_HPP

#include "internal/common.hpp"

#include <QObject>
#include <QColor>

namespace modfem {
namespace heat {

class MODFEM_HEAT_API VertexColorMapper:
	public QObject
{
		Q_OBJECT

	public:
		static const QColor INITIAL_DEFAULT_COLOR;

		Q_PROPERTY(QColor defaultColor READ defaultColor WRITE setDefaultColor NOTIFY defaultColorChanged)

		Q_PROPERTY(QVariantList map READ map WRITE setMap NOTIFY mapChanged)

		Q_PROPERTY(QByteArray colorIndices READ colorIndices WRITE setColorIndices NOTIFY colorIndicesChanged)

		Q_PROPERTY(QByteArray colorVertices READ colorVertices NOTIFY colorVerticesChanged)

		Q_PROPERTY(int count READ count NOTIFY countChanged)

		VertexColorMapper(QObject * parent = nullptr);

		QColor defaultColor() const;

		void setDefaultColor(QColor defaultColor);

		QVariantList map() const;

		void setMap(QVariantList map);

		QByteArray colorIndices() const;

		void setColorIndices(const QByteArray & colorIndices);

		QByteArray colorVertices() const;

		int count() const;

	signals:
		void defaultColorChanged();

		void mapChanged();

		void colorIndicesChanged();

		void countChanged();

		void colorVerticesChanged();

	private:
		static constexpr std::size_t COLOR_INDEX_SIZE = sizeof(int);

		static constexpr std::size_t COLOR_VECTOR_DIM = 3;

		static constexpr std::size_t COLOR_VECTOR_SIZE = sizeof(float) * COLOR_VECTOR_DIM;

		void setCount(int count);

		void updateColorVertices();

		struct Members {
			QColor defaultColor;
			QVariantList map;
			QByteArray colorIndices;
			QByteArray colorVertices;
			int count;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // VERTEXCOLORMAPPER_HPP
