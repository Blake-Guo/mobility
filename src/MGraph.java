import road.*;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.io.*;

public class MGraph {

	public RoadNetwork rnetwork;
	Map<Integer, Integer>[] mgraph;

	public MGraph() {

	}
	
	public void initMGraph(String nodefile, String edgefile)throws IOException{
		
		rnetwork = new RoadNetwork();
		rnetwork.initMap(nodefile, edgefile);
		
		mgraph = (HashMap<Integer, Integer>[]) (new HashMap<?, ?>[rnetwork.edges.size()]);
		for(int i =0;i<mgraph.length;i++)
			mgraph[i] = new HashMap<Integer, Integer>();
	}

	
	
	/**
	 * Write the generated mobility matrix to a file
	 * @param filePath
	 */
	public void graphWriter(String filePath){
		
	}
	
	/**
	 * Write the generated mobility matrix to a file in sparse form.
	 * @param filePath
	 */
	public void graphWriter_sparse(String filePath){
		
	}
	
	
	
	/**
	 * Given a trip with the start edge id and the end edge id, we add this trip to the mobility graph
	 * @param sedgeId
	 * @param eedgeId
	 */
	public void addTripToMGraph(int sedgeId, int eedgeId){
		int preval = 0;
		if(mgraph[sedgeId].containsKey(eedgeId)){
			preval = mgraph[sedgeId].get(eedgeId);
		}
		preval+=1;//add the new trip
		mgraph[sedgeId].put(eedgeId, preval);
	}

	
	/**
	 * The key function which generates the mobility feature vectors for all the edges by querying the trips in database
	 * @throws SQLException
	 */
	public void genFeatureVectors_DB() throws SQLException {
		CoordinateConversion cconvert = new CoordinateConversion();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		ResultSet rs = null;

		int day = 1;

		while (day < 3) {

			String sql = "select start_lat,start_lon,end_lat,end_lon from nytaxi "
					+ "where extract(month from start_time)=9 and "
					+ "extract(day from start_time)=" + day + ";";
			try {
				rs = statement.executeQuery(sql);
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			while(rs.next()){
				double start_lat = rs.getDouble(1);
				double start_lon = rs.getDouble(2);
				double end_lat = rs.getDouble(3);
				double end_lon = rs.getDouble(4);
				
				Point spoint = cconvert.convertLatLonToUTM(start_lat, start_lon);
				Point epoint = cconvert.convertLatLonToUTM(end_lat, end_lon);
				
				int sedgeId = rnetwork.projectPointToEdge(spoint);
				int eedgeId = rnetwork.projectPointToEdge(epoint);
				addTripToMGraph(sedgeId, eedgeId);
			}
			
			
			System.out.println("Finish processing the data in day " + day);
			day++;
		}

	}
	
	public static void main(String[] args)throws IOException, SQLException{
		MGraph mgraph = new MGraph();
		mgraph.initMGraph("data/nodes.csv", "data/edges.csv");
		mgraph.genFeatureVectors_DB();
		
	}

}
