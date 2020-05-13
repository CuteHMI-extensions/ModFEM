#ifndef H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_PROBLEM_HPP
#define H_EXTENSIONS_MODFEM_HEAT_0_INCLUDE_MODFEM_HEAT_PROBLEM_HPP

#include "internal/common.hpp"
#include "ElementData.hpp"
#include "ScalarFieldNodes.hpp"

#include <QObject>
#include <QUrl>
#include <Qt3DRender/QBuffer>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API Problem:
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

		Q_PROPERTY(ElementData * elementData READ elementData NOTIFY elementDataChanged)

		Problem(QObject * parent = nullptr);

		~Problem() override;

		QString directory() const;

		void setDirectory(const QString & directory);

		int problemId() const;

		int meshId() const;

		int fieldId() const;

		int solutionCount() const;

		int equationCount() const;

		ElementData * elementData() const;

		Q_INVOKABLE void setDirectoryFromURL(const QUrl & url);

		Q_INVOKABLE void init();

		Q_INVOKABLE void solve();

		Q_INVOKABLE void integrate();

		Q_INVOKABLE void writeParaview();

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

		void elementDataChanged();

	private:
		struct Members {
			QString directory;
			FILE * interactiveInput;
			FILE * interactiveOutput;
			int problemId;
			int meshId;
			int fieldId;
			int solutionCount;
			int equationCount;
			ElementData * elementData;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // CONTROLLER_HPP
