#include "QMLPlugin.hpp"

#include <modfem/nssupgheat/internal/common.hpp>

#include <modfem/nssupgheat/Problem.hpp>
#include <modfem/nssupgheat/ElementData.hpp>
#include <modfem/nssupgheat/AbstractColorMapper.hpp>
#include <modfem/nssupgheat/PaletteColorMapper.hpp>
#include <modfem/nssupgheat/HueColorMapper.hpp>
#include <modfem/nssupgheat/ScalarFieldNodes.hpp>
#include <modfem/nssupgheat/AbstractProbe.hpp>
#include <modfem/nssupgheat/ScalarProbe.hpp>
#include <modfem/nssupgheat/Vector3Probe.hpp>

#include <QtQml>

namespace modfem {
namespace nssupgheat {
namespace internal {

void QMLPlugin::registerTypes(const char * uri)
{
	Q_ASSERT(uri == QLatin1String("ModFEM.NSSUPGHeat"));

	qmlRegisterType<modfem::nssupgheat::Problem>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "Problem");
	qmlRegisterType<modfem::nssupgheat::PaletteColorMapper>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "PaletteColorMapper");
	qmlRegisterUncreatableType<modfem::nssupgheat::AbstractColorMapper>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "AbstractColorMapper", "Class 'modfem::heat::AbstractColorMapper' can not be instantiated from QML");
	qmlRegisterType<modfem::nssupgheat::HueColorMapper>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "HueColorMapper");
	qmlRegisterType<modfem::nssupgheat::ScalarFieldNodes>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "ScalarFieldNodes");
	qmlRegisterUncreatableType<modfem::nssupgheat::ElementData>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "ElementData", "Class 'modfem::heat::ElementData' can not be instantiated from QML");
	qmlRegisterUncreatableType<modfem::nssupgheat::AbstractProbe>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "AbstractProbe", "Class 'modfem::heat::AbstractProbe' can not be instantiated from QML");
	qmlRegisterType<modfem::nssupgheat::ScalarProbe>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "ScalarProbe");
	qmlRegisterType<modfem::nssupgheat::Vector3Probe>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "Vector3Probe");
}

}
}
}

//(c)C: Copyright © 2019-2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
