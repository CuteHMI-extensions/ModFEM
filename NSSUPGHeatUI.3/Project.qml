import QtQml 2.0

import CuteHMI.Services 2.0
import CuteHMI.SharedDatabase 0.0
import CuteHMI.DataAcquisition 0.0

QtObject {
	objectName: "project"

	property Database db: Database {
		objectName: "db"

		connectionName: "modfem::nssupgheat"
		type: "QSQLITE"
		host: "localhost"
		port: 5432
		name: "modfem_nssupgheat"
		user: "postgres"
		password: "postgres"

		/**
		  Set PostgreSQL default settings.
		  */
		function setPostgresDefaults() {
			type = "QPSQL"
			host = "localhost"
			port = 5432
			name = "postgres"
			user = "postgres"
			passowrd = "postgres"
		}

		/**
		  Set SQLite default settings.
		  */
		function setSQLiteDefaults() {
			type = "QSQLITE"
			name = "modfem_nssupgheat"
		}
	}

	property Schema schema: Schema {
		objectName: "schema"
		name: "modfem_nssupgheat"

		connectionName: "modfem::nssupgheat"
	}

	property Service service: Service {
		objectName: "service"

		name: "Database Service"
		serviceable: db
	}

}
