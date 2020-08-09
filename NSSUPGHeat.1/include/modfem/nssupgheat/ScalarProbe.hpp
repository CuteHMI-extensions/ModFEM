#ifndef SCALARPROBE_HPP
#define SCALARPROBE_HPP

#include "AbstractProbe.hpp"

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API ScalarProbe:
	public AbstractProbe

{
		Q_OBJECT

	public:
		Q_PROPERTY(double value READ value WRITE setValue NOTIFY valueChanged)

		ScalarProbe(QObject * parent = nullptr);

		double value() const;

		void setValue(double value);

	signals:
		void valueChanged();

	private:
		struct Members {
			double value;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // SCALARPROBE_HPP
