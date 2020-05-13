#include "QMLPlugin.hpp"

#include <modfem/nssupgheat/internal/common.hpp>

#include <modfem/nssupgheat/Problem.hpp>
#include <modfem/nssupgheat/ElementData.hpp>
#include <modfem/nssupgheat/VertexColorMapper.hpp>
#include <modfem/nssupgheat/ScalarFieldNodes.hpp>

#include <QtQml>

namespace modfem {
namespace nssupgheat {
namespace internal {

void QMLPlugin::registerTypes(const char * uri)
{
	Q_ASSERT(uri == QLatin1String("ModFEM.NSSUPGHeat"));

	qmlRegisterType<modfem::nssupgheat::Problem>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "Problem");
	qmlRegisterType<modfem::nssupgheat::VertexColorMapper>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "VertexColorMapper");
	qmlRegisterType<modfem::nssupgheat::ScalarFieldNodes>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "ScalarFieldNodes");
	qmlRegisterUncreatableType<modfem::nssupgheat::ElementData>(uri, MODFEM_NSSUPGHEAT_MAJOR, 0, "ElementData", "Class 'modfem::heat::ElementData' can not be instantiated from QML");
}

}
}
}

//(c)C: Copyright © 2019-2020, Michał Policht <michal@policht.pl>. All rights reserved.
//(c)C: This file is a part of CuteHMI.
//(c)C: CuteHMI is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//(c)C: CuteHMI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//(c)C: You should have received a copy of the GNU Lesser General Public License along with CuteHMI.  If not, see <https://www.gnu.org/licenses/>.
