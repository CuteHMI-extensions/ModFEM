#ifndef INTEGRATIONTHREAD_HPP
#define INTEGRATIONTHREAD_HPP

#include "internal/common.hpp"

#include <QThread>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API IntegrationThread:
	public QThread
{
		Q_OBJECT

	public:
		IntegrationThread(FILE * interactiveInput, FILE * interactiveOutput, const QString & workingDirectory, QObject * parent = nullptr);

		void run() override;

	private:
		void integrate();

		void nsSupgHeatIntegrate(char * Work_dir, FILE * Interactive_input, FILE * Interactive_output);

		void nsSupgHeatSolDiffNorm(int Current, int Old, double * sol_diff_norm_ns_supg_p, double * sol_diff_norm_heat_p);

		struct Members {
			FILE * interactiveInput;
			FILE * interactiveOutput;
			QString workingDirectory;
		};

		cutehmi::MPtr<Members> m;
};

}
}

#endif // INTEGRATIONTHREAD_HPP
