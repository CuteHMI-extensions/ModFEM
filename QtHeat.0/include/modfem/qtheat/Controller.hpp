#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_CONTROLLER_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_CONTROLLER_HPP

#include "internal/common.hpp"

#include <QObject>
#include <QUrl>

namespace modfem {
namespace qtheat {

class MODFEM_QTHEAT_API Controller:
	public QObject
{
		Q_OBJECT

	public:
		Q_PROPERTY(QString problemDirectory READ problemDirectory WRITE setProblemDirectory NOTIFY problemDirectoryChanged)

		Q_PROPERTY(int problemId READ problemId WRITE setProblemId NOTIFY problemIdChanged)

		Q_PROPERTY(int meshId READ meshId WRITE setMeshId NOTIFY meshIdChanged)

		Q_PROPERTY(int fieldId READ fieldId WRITE setFieldId NOTIFY fieldIdChanged)

		Q_PROPERTY(int solutionCount READ solutionCount WRITE setSolutionCount NOTIFY solutionCountChanged)

		Q_PROPERTY(int equationCount READ equationCount WRITE setEquationCount NOTIFY equationCountChanged)

		Controller(QObject * parent = nullptr);

		~Controller() override;

		QString problemDirectory() const;

		void setProblemDirectory(const QString & problemDirectory);

		int problemId() const;

		int meshId() const;

		int fieldId() const;

		int solutionCount() const;

		int equationCount() const;

		Q_INVOKABLE void setProblemDirectoryFromURL(const QUrl & url);

		Q_INVOKABLE void init();

	protected slots:
		void setProblemId(int problemId);

		void setMeshId(int meshId);

		void setFieldId(int fieldId);

		void setSolutionCount(int solutionCount);

		void setEquationCount(int equationCount);

	signals:
		void problemDirectoryChanged();

		void problemIdChanged();

		void meshIdChanged();

		void fieldIdChanged();

		void solutionCountChanged();

		void equationCountChanged();

	private:
		struct Members {
			QString problemDirectory;
			int problemId;
			int meshId;
			int fieldId;
			int solutionCount;
			int equationCount;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // CONTROLLER_HPP
