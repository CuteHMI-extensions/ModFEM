#ifndef VECTOR3PROBE_HPP
#define VECTOR3PROBE_HPP

#include "AbstractProbe.hpp"

#include <QVector3D>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API Vector3Probe:
	public AbstractProbe
{
		Q_OBJECT

	public:
		Q_PROPERTY(QVector3D value READ value WRITE setValue NOTIFY valueChanged)

		Vector3Probe(QObject * parent = nullptr);

		QVector3D value() const;

		void setValue(const QVector3D & value);

	signals:
		void valueChanged();

	private:
		struct Members {
			QVector3D value;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // VECTORPROBE_HPP
