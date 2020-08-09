#ifndef ABSTRACTPROBE_HPP
#define ABSTRACTPROBE_HPP

#include "internal/common.hpp"

#include <cutehmi/MPtr.hpp>

#include <QObject>
#include <QVector3D>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API AbstractProbe:
	public QObject
{
		Q_OBJECT

	public:
		Q_PROPERTY(QVector3D position READ position WRITE setPosition NOTIFY positionChanged)

		Q_PROPERTY(QString field READ field WRITE setField NOTIFY fieldChanged)

		AbstractProbe(QObject * parent = nullptr);

		QVector3D position() const;

		void setPosition(const QVector3D position);

		QString field() const;

		void setField(const QString & field);

	signals:
		void positionChanged();

		void fieldChanged();

	private:
		struct Members {
			QVector3D position;
			QString field;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // ABSTRACTPROBE_HPP
