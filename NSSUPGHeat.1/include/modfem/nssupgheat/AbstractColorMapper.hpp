#ifndef ABSTRACTCOLORMAPPER_HPP
#define ABSTRACTCOLORMAPPER_HPP

#include "internal/common.hpp"

#include <QObject>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API AbstractColorMapper:
	public QObject
{
		Q_OBJECT

	public:
		enum InputType {
			INPUT_FLOAT,
			INPUT_DOUBLE,
			INPUT_INT
		};
		Q_ENUM(InputType)

		static constexpr InputType INITIAL_INPUT_TYPE = INPUT_FLOAT;

		Q_PROPERTY(InputType inputType READ inputType WRITE setInputType NOTIFY inputTypeChanged)

		Q_PROPERTY(QByteArray input READ input WRITE setInput NOTIFY inputChanged)

		Q_PROPERTY(QByteArray output READ output NOTIFY outputChanged)

		Q_PROPERTY(int count READ count NOTIFY countChanged)

		InputType inputType() const;

		void setInputType(InputType type);

		QByteArray input() const;

		void setInput(const QByteArray & input);

		QByteArray output() const;

		int count() const;

	signals:
		void inputTypeChanged();

		void inputChanged();

		void outputChanged();

		void countChanged();

	protected:
		static constexpr std::size_t OUTPUT_DIM = 3;

		static constexpr std::size_t OUTPUT_SIZE = sizeof(float) * OUTPUT_DIM;

		AbstractColorMapper(QObject * parent = nullptr);

		void setCount(int count);

		std::size_t inputTypeSize() const;

		QByteArray & input();

		QByteArray & output();

	private:
		struct Members {
			InputType inputType;
			QByteArray input;
			QByteArray output;
			int count;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // ABSTRACTCOLORMAPPER_HPP
