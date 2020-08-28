#ifndef INTEGRATIONDATA_HPP
#define INTEGRATIONDATA_HPP

#include "internal/common.hpp"

#include <QObject>
#include <QReadWriteLock>

namespace modfem {
namespace nssupgheat {

class MODFEM_NSSUPGHEAT_API IntegrationData:
	public QObject
{
		Q_OBJECT

	public:
		static constexpr bool INITIAL_REAL_TIME_SIMULATION = true;

		static constexpr double INITIAL_REQUESTED_TIME_STEP = 0.0;

		Q_PROPERTY(bool realTimeSimulation READ realTimeSimulation WRITE setRealTimeSimulation NOTIFY realTimeSimulationChanged)

		Q_PROPERTY(double requestedTimeStep READ requestedTimeStep WRITE setRequestedTimeStep NOTIFY requestedTimeStepChanged)

		Q_PROPERTY(double currentTimeStep READ currentTimeStep NOTIFY currentTimeStepChanged)

		Q_PROPERTY(double realTime READ realTime NOTIFY realTimeChanged)

		Q_PROPERTY(double simulationTime READ simulationTime NOTIFY simulationTimeChanged)

		IntegrationData(QObject * parent = nullptr);

		bool realTimeSimulation() const;

		void setRealTimeSimulation(bool realTimeSimulation);

		double requestedTimeStep() const;

		void setRequestedTimeStep(double requestedTimeStep);

		double currentTimeStep() const;

		double realTime() const;

		double simulationTime() const;

	public slots:
		void setCurrentTimeStep(double currentTimeStep);

		void setSimulationTime(double simulationTime);

		void setRealTime(double realTime);

	signals:
		void realTimeSimulationChanged();

		void requestedTimeStepChanged();

		void currentTimeStepChanged();

		void realTimeChanged();

		void simulationTimeChanged();

	private:
		struct Members {
			mutable QReadWriteLock lock;
			bool realTimeSimulation;
			double requestedTimeStep;
			double currentTimeStep;
			double realTime;
			double simulationTime;

			Members():
				realTimeSimulation(INITIAL_REAL_TIME_SIMULATION),
				requestedTimeStep(INITIAL_REQUESTED_TIME_STEP),
				currentTimeStep(std::numeric_limits<double>::quiet_NaN()),
				realTime(std::numeric_limits<double>::quiet_NaN()),
				simulationTime(std::numeric_limits<double>::quiet_NaN())
			{
			}
		};

		cutehmi::MPtr<Members> m;

};

}
}

#endif // INTEGRATIONDATA_HPP
