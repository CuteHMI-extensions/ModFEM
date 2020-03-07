#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_FACEDATA_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_FACEDATA_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace qtheat {

class FaceData:
	public QObject
{
		Q_OBJECT

	public:
		explicit FaceData(QObject * parent = nullptr);

		Q_PROPERTY(int triangleCount READ triangleCount NOTIFY triangleCountChanged)

		Q_PROPERTY(QByteArray triangleIndices READ triangleIndices NOTIFY triangleIndicesChanged)

		Q_PROPERTY(QByteArray triangleCoords READ triangleCoords NOTIFY triangleCoordsChanged)

		Q_PROPERTY(QByteArray triangleNormals READ triangleNormals NOTIFY triangleNormalsChanged)

		Q_PROPERTY(int quadCount READ quadCount NOTIFY quadCountChanged)

		Q_PROPERTY(QByteArray quadIndices READ quadIndices NOTIFY quadIndicesChanged)

		Q_PROPERTY(QByteArray quadCoords READ quadCoords NOTIFY quadCoordsChanged)

		int triangleCount() const;

		QByteArray triangleIndices() const;

		QByteArray triangleCoords() const;

		QByteArray triangleNormals() const;

		int quadCount() const;

		QByteArray quadIndices() const;

		QByteArray quadCoords();

	public slots:
		void init(int meshId);

	signals:
		void triangleCountChanged();

		void triangleIndicesChanged();

		void triangleCoordsChanged();

		void triangleNormalsChanged();

		void quadCountChanged();

		void quadIndicesChanged();

		void quadCoordsChanged();

	protected slots:
		void setTriangleCount(int count);

		void setQuadCount(int count);

	private:
		struct Members {
			int triangleCount;
			QByteArray triangleIndices;
			QByteArray triangleCoords;
			QByteArray triangleNormals;
			int quadCount;
			QByteArray quadIndices;
			QByteArray quadCoords;
		};

		cutehmi::MPtr<Members> m;

};

}
}

#endif // FACEDATA_HPP
