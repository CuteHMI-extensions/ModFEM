#ifndef H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_PROBLEM_HPP
#define H_EXTENSIONS_MODFEM_QTHEAT_0_INCLUDE_MODFEM_QTHEAT_PROBLEM_HPP

#include "internal/common.hpp"
#include "Mesh.hpp"

#include <QObject>
#include <QUrl>
#include <Qt3DRender/QBuffer>

namespace modfem {
namespace qtheat {

class MODFEM_QTHEAT_API Problem:
	public QObject
{
		Q_OBJECT

	public:
		Q_PROPERTY(QString directory READ directory WRITE setDirectory NOTIFY directoryChanged)

		Q_PROPERTY(int problemId READ problemId WRITE setProblemId NOTIFY problemIdChanged)

		Q_PROPERTY(int meshId READ meshId WRITE setMeshId NOTIFY meshIdChanged)

		Q_PROPERTY(int fieldId READ fieldId WRITE setFieldId NOTIFY fieldIdChanged)

		Q_PROPERTY(int solutionCount READ solutionCount WRITE setSolutionCount NOTIFY solutionCountChanged)

		Q_PROPERTY(int equationCount READ equationCount WRITE setEquationCount NOTIFY equationCountChanged)

		Q_PROPERTY(Mesh * mesh READ mesh NOTIFY meshChanged)

		Problem(QObject * parent = nullptr);

		~Problem() override;

		QString directory() const;

		void setDirectory(const QString & directory);

		int problemId() const;

		int meshId() const;

		int fieldId() const;

		int solutionCount() const;

		int equationCount() const;

		Mesh * mesh() const;

		Q_INVOKABLE void setDirectoryFromURL(const QUrl & url);

		Q_INVOKABLE void init();

	protected slots:
		void setProblemId(int problemId);

		void setMeshId(int meshId);

		void setFieldId(int fieldId);

		void setSolutionCount(int solutionCount);

		void setEquationCount(int equationCount);

	signals:
		void directoryChanged();

		void problemIdChanged();

		void meshIdChanged();

		void fieldIdChanged();

		void solutionCountChanged();

		void equationCountChanged();

		void meshChanged();

	private:
		struct Members {
			QString directory;
			int problemId;
			int meshId;
			int fieldId;
			int solutionCount;
			int equationCount;
			Mesh * mesh;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // CONTROLLER_HPP
