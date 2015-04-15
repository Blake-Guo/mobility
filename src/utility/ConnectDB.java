package utility;


import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
 

public class ConnectDB {
	public static Connection ConnectDatabase(String dbname) {
		try {
			Class.forName("org.postgresql.Driver");
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}

		System.out.println("PostgreSQL JDBC Driver Registered!");

		Connection connection = null;

		try {
			connection = DriverManager.getConnection(
					"jdbc:postgresql://107.170.82.229:5432/" + dbname,
					"postg", "abcdefg");
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}

		System.out.println("Successfully connected to the data base");
		return connection;
	}
}
