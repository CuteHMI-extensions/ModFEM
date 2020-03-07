#include "QMLPlugin.hpp"

#include <modfem/qtheat/Problem.hpp>
#include <modfem/qtheat/Mesh.hpp>
#include <modfem/qtheat/FaceData.hpp>

#include <QtQml>

namespace modfem {
namespace qtheat {
namespace internal {

void QMLPlugin::registerTypes(const char * uri)
{
	Q_ASSERT(uri == QLatin1String("ModFEM.QtHeat"));

	qmlRegisterType<modfem::qtheat::Problem>(uri, MODFEM_QTHEAT_MAJOR, 0, "Problem");
	qmlRegisterType<modfem::qtheat::Mesh>(uri, MODFEM_QTHEAT_MAJOR, 0, "Mesh");
	qmlRegisterUncreatableType<modfem::qtheat::FaceData>(uri, MODFEM_QTHEAT_MAJOR, 0, "FaceData", "Class 'modfem::qtheat::FaceData' can not be instantiated from QML");
}

}
}
}

//(c)C: Copyright © 2019-2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
