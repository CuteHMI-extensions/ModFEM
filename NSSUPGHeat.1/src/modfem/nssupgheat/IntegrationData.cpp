#include "../../../include/modfem/nssupgheat/IntegrationData.hpp"

namespace modfem {
namespace nssupgheat {

constexpr bool IntegrationData::INITIAL_REAL_TIME_SIMULATION;
constexpr double IntegrationData::INITIAL_REQUESTED_TIME_STEP;

IntegrationData::IntegrationData(QObject * parent):
	QObject(parent),
	m(new Members)
{
}

bool IntegrationData::realTimeSimulation() const
{
	QReadLocker locker(& m->lock);

	return m->realTimeSimulation;
}

void IntegrationData::setRealTimeSimulation(bool realTime)
{
	QReadLocker locker(& m->lock);

	if (m->realTimeSimulation != realTime) {
		m->lock.unlock();
		m->lock.lockForWrite();
		m->realTimeSimulation = realTime;
		m->lock.unlock();
		emit realTimeSimulationChanged();
	}
}

double IntegrationData::requestedTimeStep() const
{
	QReadLocker locker(& m->lock);

	return m->requestedTimeStep;
}

void IntegrationData::setRequestedTimeStep(double timeStep)
{
	QReadLocker locker(& m->lock);

	if (m->requestedTimeStep != timeStep) {
		m->lock.unlock();
		m->lock.lockForWrite();
		m->requestedTimeStep = timeStep;
		m->lock.unlock();
		emit requestedTimeStepChanged();
	}
}

double IntegrationData::currentTimeStep() const
{
	return m->currentTimeStep;
}

double IntegrationData::realTime() const
{
	return m->realTime;
}

double IntegrationData::simulationTime() const
{
	return m->simulationTime;
}

void IntegrationData::setCurrentTimeStep(double currentTimeStep)
{
	if (m->currentTimeStep != currentTimeStep) {
		m->currentTimeStep = currentTimeStep;
		emit currentTimeStepChanged();
	}
}

void IntegrationData::setSimulationTime(double simulationTime)
{
	if (m->simulationTime != simulationTime) {
		m->simulationTime = simulationTime;
		emit simulationTimeChanged();
	}
}

void IntegrationData::setRealTime(double realTime)
{
	if (m->realTime != realTime) {
		m->realTime = realTime;
		emit realTimeChanged();
	}
}

}
}
