import road.*;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.*;
import java.io.*;
import java.text.*;

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
	public void writeMGraph_dense(String filePath){
		
	}
	
	/**
	 * Write the generated mobility matrix to a file in sparse form.
	 * @param filePath
	 */
	public void writeMGraph_sparse(String filePath)throws IOException{
		BufferedWriter bwriter = new BufferedWriter(new FileWriter(filePath));
		for(int i=0;i<rnetwork.edges.size();i++){
			StringBuilder record = new StringBuilder();

			//bwriter.write(i);
			record.append(i);
			
			for(Map.Entry<Integer, Integer> entry : mgraph[i].entrySet()){
				//bwriter.write("\t" + entry.getKey()+":" + entry.getValue());
				record.append('\t');
				record.append(entry.getKey());
				record.append(':');
				record.append(entry.getValue());
			}
			//record.append('\n');
			bwriter.write(record.toString());
			bwriter.newLine();
		}
		bwriter.close();
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
	 * @throws ParseException 
	 */
	public void genFeatureVectors_DB() throws SQLException, ParseException {
		CoordinateConversion cconvert = new CoordinateConversion();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		ResultSet rs = null;

		//int day = 1;
		int incorrectTrip = 0;
		SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
		String bTimeStr = new String("2013-09-01");
		String eTimeStr = new String("2013-10-31");
		
		Calendar c_endTime = Calendar.getInstance();
		c_endTime.setTime(dateFormat.parse(eTimeStr));
		
		Calendar c_time = Calendar.getInstance();
		c_time.setTime(dateFormat.parse(bTimeStr));

		while (c_time.after(c_endTime) == false) {

			String sql = "select start_lat,start_lon,end_lat,end_lon from nytaxi "
					+ "where extract(month from start_time)=" + (c_time.get(Calendar.MONTH)+1) +
					" and "
					+ "extract(day from start_time)=" + c_time.get(Calendar.DAY_OF_MONTH)  + ";";
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
				if(sedgeId == -1 || eedgeId == -1)
				{
					incorrectTrip++;
					continue;
				}
				addTripToMGraph(sedgeId, eedgeId);
			}
			rs.close();
			
			System.out.println("Finish processing the data in day " + dateFormat.format(c_time.getTime()));
			c_time.add(Calendar.DATE, 1);
		}
		
		System.out.println("The number of incorrect trip is: " + incorrectTrip);

	}
	
	public static void main(String[] args)throws IOException, SQLException, ParseException{
		MGraph mgraph = new MGraph();
		mgraph.initMGraph("data/nodes.csv", "data/edges.csv");
		mgraph.genFeatureVectors_DB();
		mgraph.writeMGraph_sparse("sparse_mgraph.csv");
	}

}
