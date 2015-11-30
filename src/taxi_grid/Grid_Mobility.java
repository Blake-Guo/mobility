package taxi_grid;

import utility.ConnectDB;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TimeZone;

import Geo.Point_latlon;
import Geo.Polygon;
import Geo.Point;
import road.CoordinateConversion;

/*
 * Here we assume the boundary is a parallelogram.
 * The boundary is decided by four parameters: (1) four boundPoints, starting from south_west counterclockwise; (2) the width and height of the whole boundary 
 */
public class Grid_Mobility {
	
	Point[] boundPoints;
	Point_latlon[] boundGeoPoints;
	double grid_width;// this is the width of a grid, and the unit metric is kilometer
	double grid_height;// this is the height of a grid, and the unit metric kilometer
	
	double boundary_width;// the width of the boundary
	double boundary_height;// the height of the boundary
	Polygon[] grid_boundaries;
	public int grid_xnumber;
	public int grid_ynumber;
	int grid_number;
	int[][][] mtensor;//the mode of mtensor is reverse of the real mobility tensor. Here mode 1: time, mode 2: original, mode 3: destination
	final int tensor_time_period = 24;
	int ncluster; // the number of cluster
	Map<Integer, Integer> grid2cluster;//store the cluster id of each grid
	Map<Integer, ArrayList<Integer>> clusteredGrids;//store the grid ids of each cluster
	
	public Grid_Mobility(){
		//Step0: The boundary of Manhatten, (-74.029999,40.703546), (-73.974380,40.703546), (-73.924942,40.791459),(-73.98056,40.791459)
		CoordinateConversion convert = new CoordinateConversion();
		boundGeoPoints = new Point_latlon[4];
		boundGeoPoints[0] = new Point_latlon(40.703546, -74.029999);
		boundGeoPoints[1] = new Point_latlon(40.703546, -73.974380);
		boundGeoPoints[2] = new Point_latlon(40.791459, -73.924942);
		boundGeoPoints[3] = new Point_latlon(40.791459, -73.98056);
		
		boundPoints = new Point[4];
		boundPoints[0] = convert.convertLatLonToUTM(40.703546, -74.029999);
		boundPoints[1] = convert.convertLatLonToUTM(40.703546, -73.974380);
		boundPoints[2] = convert.convertLatLonToUTM(40.791459, -73.924942);
		boundPoints[3] = convert.convertLatLonToUTM(40.791459, -73.98056);
		
		//Step1: initialize the grids
		grid_width = 800;//meter
		grid_height = 800;//meter
		boundary_width = boundPoints[1].x - boundPoints[0].x;
		boundary_height = boundPoints[2].y - boundPoints[1].y;
		grid_xnumber = (int)(boundary_width/grid_width);
		grid_ynumber = (int)(boundary_height/grid_height);
		grid_number = grid_xnumber * grid_ynumber;
		
	}
	
	
	
	public Grid_Mobility(String clusterFile) throws IOException{
		//Step0: The boundary of Manhatten, (-74.029999,40.703546), (-73.974380,40.703546), (-73.924942,40.791459),(-73.98056,40.791459)
		CoordinateConversion convert = new CoordinateConversion();
		boundGeoPoints = new Point_latlon[4];
		boundGeoPoints[0] = new Point_latlon(40.703546, -74.029999);
		boundGeoPoints[1] = new Point_latlon(40.703546, -73.974380);
		boundGeoPoints[2] = new Point_latlon(40.791459, -73.924942);
		boundGeoPoints[3] = new Point_latlon(40.791459, -73.98056);
		
		boundPoints = new Point[4];
		boundPoints[0] = convert.convertLatLonToUTM(40.703546, -74.029999);
		boundPoints[1] = convert.convertLatLonToUTM(40.703546, -73.974380);
		boundPoints[2] = convert.convertLatLonToUTM(40.791459, -73.924942);
		boundPoints[3] = convert.convertLatLonToUTM(40.791459, -73.98056);
		
		//Step1: initialize the grids
		grid_width = 800;//meter
		grid_height = 800;//meter
		
		boundary_width = boundPoints[1].x - boundPoints[0].x;
		boundary_height = boundPoints[2].y - boundPoints[1].y;
		grid_xnumber = (int)(boundary_width/grid_width);
		grid_ynumber = (int)(boundary_height/grid_height);
		grid_number = grid_xnumber * grid_ynumber;
		
		grid_boundaries = new Polygon[grid_xnumber * grid_ynumber];
		int grid_index = 0;
		double x_level_offset = (boundPoints[2].x - boundPoints[1].x) / grid_ynumber;
		for (int i = 0; i < grid_ynumber; i++)
			for (int j = 0; j < grid_xnumber; j++) {
				
				grid_boundaries[grid_index] = new Polygon();
				
				Point lbpoint = new Point(boundPoints[0].x + i
						* x_level_offset + j * grid_width, boundPoints[0].y + i
						* grid_height);
				grid_boundaries[grid_index].addPoint(lbpoint);
				
				Point rbpoint = new Point(lbpoint.x + grid_width, lbpoint.y);
				grid_boundaries[grid_index].addPoint(rbpoint);
				
				Point rupoint = new Point(lbpoint.x + grid_width+ x_level_offset,
						lbpoint.y + grid_height);
				grid_boundaries[grid_index].addPoint(rupoint);
				
				Point lupoint = new Point(rupoint.x - grid_width, rupoint.y);
				grid_boundaries[grid_index].addPoint(lupoint);
				
				grid_index++;
			}
		
		
		//Step2: Read in the cluster file and build the grid to cluster
		BufferedReader br = new BufferedReader(new FileReader(clusterFile));
		ncluster = 0;
		Set<Integer> clusterIds = new HashSet<Integer>();
		grid2cluster = new HashMap<Integer,Integer>();
		clusteredGrids = new HashMap<Integer, ArrayList<Integer>>();
		String line = "";
		
		while((line = br.readLine()) != null){
			String[] strs = line.split(" ");
			int gridId = Integer.valueOf(strs[0]);
			int cId = Integer.valueOf(strs[1]);
			
			if(clusterIds.contains(cId) == false){
				ncluster++;
				clusterIds.add(cId);
			}
			
			grid2cluster.put(gridId, cId);
		}
		
		
		//step3: compute the grids of each cluster
		for(Map.Entry<Integer, Integer> entry : grid2cluster.entrySet()){
			if(clusteredGrids.containsKey(entry.getValue()) == false){
				int cid = entry.getValue();//cluster id
				clusteredGrids.put(cid, new ArrayList<Integer>());
			}
			clusteredGrids.get(entry.getValue()).add(entry.getKey());
		}
		
		System.out.println("Finished Constructing");
	}
	
	
	
	
	public void outputGridsToFile(String ofile, int[] cids) throws IOException{
		BufferedWriter br = new BufferedWriter(new FileWriter(ofile));
	
		br.write("cInd\t"+"cId\t" + "gId\t" + "latitude\t"+ "longitude\n");
		int cIndex = 1;
		for(int i=0;i<cids.length;i++) // for each cluster
		{
			//step0: get the grids within the specific cluster
			ArrayList<Integer> grids = clusteredGrids.get(cids[i]);
			
			for(int j=0;j<grids.size();j++){// for each grid
				ArrayList<Point_latlon> point_latlons = getGridGeoBoundary(grids.get(j));
				
				for(Point_latlon p : point_latlons){
					String str = Integer.toString(cIndex) + "\t" + Integer.toString(cids[i]) + "\t" + Integer.toString(grids.get(j)) + "\t";
					str += Double.toString(p.lat) + "\t";
					str += Double.toString(p.lon) + "\n";
					br.write(str);
				}
			}
			cIndex++;
		}
		br.close();
		System.out.println("Finish writing cluster geo boundary");
	}
	
	
	/**
	 * Convert the given polygon to a WKT format(String).
	 * @param poly
	 * @return
	 */
	public String toPostGISWKT(Polygon poly){
		String strbeg = "ST_GeomFromText('POLYGON((";
		String strend = "))',92030)";
		String strmid = "";
		for(Point p : poly.points){
			strmid += p.x + " ";
			strmid += p.y + ",";
		}
		strmid += poly.points.get(0).x + " ";
		strmid += poly.points.get(0).y;
		
		String str = strbeg + strmid + strend;
		
		return str;
	}
	
	
	/**
	 * Get the four geo points of a specific grid
	 * @param gIndex
	 * @return
	 */
	public ArrayList<Point_latlon> getGridGeoBoundary(int gIndex){
		ArrayList<Point_latlon> point_latlons = new ArrayList<Point_latlon>();
		double grid_lon = (boundGeoPoints[1].lon - boundGeoPoints[0].lon) / grid_xnumber;
		double grid_lat = (boundGeoPoints[2].lat - boundGeoPoints[1].lat) / grid_ynumber;

		double lon_level_offset = (boundGeoPoints[2].lon - boundGeoPoints[1].lon) / grid_ynumber;
		int i = gIndex / grid_xnumber;
		int j = gIndex % grid_xnumber;

		Point_latlon lbpoint = new Point_latlon(boundGeoPoints[0].lat + i
				* grid_lat,boundGeoPoints[0].lon + i
				* lon_level_offset + j * grid_lon);
		point_latlons.add(lbpoint);
		
		Point_latlon rbpoint = new Point_latlon(lbpoint.lat, lbpoint.lon + grid_lon);
		point_latlons.add(rbpoint);

		Point_latlon rupoint = new Point_latlon(lbpoint.lat + grid_lat, rbpoint.lon + lon_level_offset);
		point_latlons.add(rupoint);
		
		Point_latlon lupoint = new Point_latlon(rupoint.lat, rupoint.lon - grid_lon);
		point_latlons.add(lupoint);
		
		
		return point_latlons;
	}
	
	

	
	/**
	 * Get average 24 hours out-flow of a given cluster
	 * @param cId
	 * @return
	 * @throws SQLException
	 */
	public ArrayList<Integer> temporalFlow(int cId) throws SQLException{
		ArrayList<Integer> flows = new ArrayList<Integer>();
		Map<String, Integer> t2v = new HashMap<String, Integer>();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		
		System.out.println("Cluster ID:" + cId);

		ArrayList<Polygon> polyset = new ArrayList<Polygon>();
		for(int i=0;i<grid_number;i++)
			if(grid2cluster.get(i) == cId)
			{
				System.out.print("grid:" + i +",");
				polyset.add(grid_boundaries[i]);
			}
		System.out.println();
		
		for(int i=0;i<polyset.size();i++){
			
			String str_grid = toPostGISWKT(polyset.get(i));
				// we ignore the time interval currently
				String sql = "select extract(hour from start_time) as hour, count(*) as number from nytaxi"  
						+ " where ST_within(start_location," + str_grid + ") "
						+ " and extract(dow from start_time) != 0 and extract(dow from start_time) != 6 and extract(dow from start_time) != 5"
						+ " group by extract(hour from start_time)"
						+ " order by extract(hour from start_time)";
	
				ResultSet rs = statement.executeQuery(sql);
						
				//Go over each trip record
				while (rs.next()) {
					String h = rs.getString(1);
					int c = rs.getInt(2);
					String key = h;
					t2v.put(key, c);
				}
				System.out.println("Finish Grid Index:" + i);
			}	
				
			System.out.println("Finish Quering Data");
			
			for(int i=0;i<24;i++)
			{
				int val = t2v.get(Integer.toString(i));
				System.out.print(val + ",");
				flows.add(val);
			}
			System.out.println();
			
			return flows;
		}
	
	/**
	 * Get the out-flow volumes of a given cluster id at a specific time interval
	 * @param cId
	 * @param hour
	 * @return
	 * @throws SQLException 
	 */
	public ArrayList<Integer> flowDistribution(int cId, int hour) throws SQLException{
		ArrayList<Integer> flows = new ArrayList<Integer>();
		Map<String, Integer> t2v = new HashMap<String, Integer>();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		
		ArrayList<Polygon> polyset = new ArrayList<Polygon>();
		for(int i=0;i<grid_number;i++)
			if(grid2cluster.get(i) == cId)
			{
				System.out.print("grid:" + i +",");
				polyset.add(grid_boundaries[i]);
			}
		System.out.println();
		
		for(int i=0;i<polyset.size();i++){
			
			String str_grid = toPostGISWKT(polyset.get(i));
				// we ignore the time interval currently
				String sql = "select extract(month from start_time) as month, extract(day from start_time) as day, count(*) as number from nytaxi"  
						+ " where ST_within(start_location," + str_grid + ") and extract(hour from start_time) =" + Integer.toString(hour)
						+ " and extract(dow from start_time) != 0 and extract(dow from start_time) != 6 and extract(dow from start_time) != 5"
						+ " group by extract(month from start_time), extract(day from start_time)"
						+ " order by extract(month from start_time), extract(day from start_time)";
	
				ResultSet rs = statement.executeQuery(sql);
						
				//Go over each trip record
				while (rs.next()) {
					String m = rs.getString(1);
					String d = rs.getString(2);
					int c = rs.getInt(3);
					String key = m + "_" + d;
					int tmp_val = 0;
					if(t2v.containsKey(key)){
						tmp_val = t2v.get(key);
					}
					tmp_val += c;
					t2v.put(key, tmp_val);
				}
				System.out.println("Finish Grid Index:" + i);
			}	
				
			System.out.println("Finish Quering Data");
			
			for(Map.Entry<String, Integer> entry : t2v.entrySet())
			{
				System.out.print(entry.getValue() + ",");
				flows.add(entry.getValue());
			}
			
			return flows;
		}
	
	
	/**
	 * Given a point, check which grid/neighborhood it is located, and return the index.
	 * @param p
	 * @return
	 */
	public int locatePointInGridNeighbor(Point p){
		for(int i=0;i<grid_number;i++)
		{
			//Polygon grid = getGridBoundaryBasedOnIndex(i);
			Polygon grid = grid_boundaries[i];
			if(grid.contain(p))
				return i;
		}
				
		return -1;
	}
	
//	public Polygon getGridBoundaryBasedOnIndex(int index){
//		Polygon grid_boundary = new Polygon();
//		double x_level_offset = (boundPoints[2].x - boundPoints[1].x) / grid_ynumber;
//		int i = index / grid_xnumber;
//		int j = index % grid_xnumber;
//
//		Point lbpoint = new Point(boundPoints[0].x + i
//				* x_level_offset + j * grid_width, boundPoints[0].y + i
//				* grid_height);
//		grid_boundary.addPoint(lbpoint);
//		
//		Point rbpoint = new Point(lbpoint.x + grid_width, lbpoint.y);
//		grid_boundary.addPoint(rbpoint);
//		
//		Point rupoint = new Point(lbpoint.x + grid_width+ x_level_offset,
//				lbpoint.y + grid_height);
//		grid_boundary.addPoint(rupoint);
//		
//		Point lupoint = new Point(rupoint.x - grid_width, rupoint.y);
//		grid_boundary.addPoint(lupoint);
//		
//		return grid_boundary;
//	}
	
	
	/**
	 * Convert the java calendar to a timestamp string for PostgreSQL, where time_zone_str is for PostgreSQL time zone.
	 * @param time
	 * @param time_zone_str
	 * @return
	 */
	public String transfer2TimeStampStr(Calendar time, String time_zone_str){
		
		SimpleDateFormat f = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		f.setTimeZone(time.getTimeZone());
		String timestamp = f.format(time.getTime());
		//timestamp = "timestamp with time zone '" + timestamp + "'";
		timestamp = "'" + timestamp + time_zone_str + "'";		
		return timestamp;
	}
	
	
	/**
	 * Convert a wkt to Point_latlon, wkt is like  "POINT(-73.981239 40.74427)"
	 * @param wkt
	 * @return
	 */
	public Point_latlon convertWKT2PointLatlon(String wkt){
		wkt = wkt.substring(6,wkt.length()-1);
		String[] xy = wkt.split(" ");
		Point_latlon p = new Point_latlon(Double.valueOf(xy[1]), Double.valueOf(xy[0]));
		return p;
	}
	
	/**
	 * Read in the latent features(probably origin, destination, time). Each row is a record and each column is a feature.
	 * @param file
	 * @param sep
	 * @return
	 * @throws IOException 
	 */
	public ArrayList<ArrayList<Float>> readInFeatures(String file, String sep) throws IOException{
		ArrayList<ArrayList<Float>> latent_features = new ArrayList<ArrayList<Float>>();
		BufferedReader bread = new BufferedReader(new FileReader(file));
		String line = "";
		while( (line = bread.readLine()) != null ){
			String[] strs = line.split(sep);
			ArrayList<Float> new_record = new ArrayList<Float>();
			for(String s : strs){
				float val = Float.valueOf(s);
				new_record.add(val);
			}
			
			latent_features.add(new_record);
		}
		bread.close();
		return latent_features;
	}
	
	
	/**
	 * Concatenate all latent features into a string.
	 * @param dest_feature
	 * @param orig_feature
	 * @param time_feature
	 * @param sep
	 * @return
	 */
	public String convertLatentFeature2String(ArrayList<Float> orig_feature,ArrayList<Float> dest_feature, ArrayList<Float> time_feature, String sep){
		String line = "";
		
		line = orig_feature.get(0).toString();
		for(int i=1;i<orig_feature.size();i++)
			line += sep + orig_feature.get(i).toString();
		
		for(int i=0;i<dest_feature.size();i++)
			line += sep + dest_feature.get(i).toString();
		
		for(int i=0;i<time_feature.size();i++)
			line += sep + time_feature.get(i).toString();
		
		return line;
	}
	
	
	/**
	 * Convert all latent features into a string.
	 * @param feature1
	 * @param feature2
	 * @param sep
	 * @return
	 */
	public String convertLatentFeature2String(ArrayList<Float> feature1, ArrayList<Float> feature2, String sep){
		String line = "";
		
		line = feature1.get(0).toString();
		for(int i=1;i<feature1.size();i++)
			line += sep + feature1.get(i).toString();
		
		for(int i=0;i<feature2.size();i++)
			line += sep + feature2.get(i).toString();
		
		return line;
	}
	
	

	public void gen_gp_inputs(String o_flow_file) throws SQLException, IOException{
		
		//Step0: Set the begining time and end time of training data.
		Calendar btime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		Calendar etime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		btime.set(2013, 8, 3, 0, 0, 0);//the month begins from zero.
		etime.set(2013, 9, 20, 0, 0, 0);//the month begins from zero.
		String sep = " ";
		
		//step1: read the latent feature vectors from several given files. 
		//ArrayList<ArrayList<Float>> dest_features = readInFeatures(i_d_file,sep);
		//ArrayList<ArrayList<Float>> orig_features = readInFeatures(i_o_file,sep);
		//ArrayList<ArrayList<Float>> time_features = readInFeatures(i_t_file,sep);
		
		int nneigh = ncluster;
		System.out.println("Number of cluster:" + ncluster);

		//Step2: Get the trip records from the database
        int[][] hour_flow = new int[nneigh][];
        for(int i=0;i<nneigh;i++)
        	hour_flow[i] = new int[nneigh];
        
		BufferedWriter o_flow_Writer = new BufferedWriter(new FileWriter(o_flow_file));
		Random rand = new Random();
        //BufferedWriter ofeature_Writer = new BufferedWriter(new FileWriter(o_index_file));
        
		//BufferedWriter dest_valWriter = new BufferedWriter(new FileWriter(o_dval_file));
        //BufferedWriter dest_featureWriter = new BufferedWriter(new FileWriter(o_dfeature_file));
        
		CoordinateConversion convert = new CoordinateConversion();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		Calendar ctime = (Calendar)btime.clone();
		while(ctime.before(etime)) {
			if(ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SUNDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SATURDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.FRIDAY){
				ctime.add(Calendar.HOUR_OF_DAY, 1);
				continue;
			}
			
			for(int i=0;i<nneigh;i++)
				for(int j=0;j<nneigh;j++)
					hour_flow[i][j] = 0;

			//String str_grid_rect = transfer2PostGISPolygon(llpoint, rupoint);
			Calendar interval_beg = (Calendar)ctime.clone();
			Calendar interval_end = (Calendar)ctime.clone();
			interval_end.add(Calendar.HOUR, 1);
			String interval_beg_str = transfer2TimeStampStr(interval_beg,"EDT");
			String interval_end_str = transfer2TimeStampStr(interval_end,"EDT");
			//System.out.println(interval_beg_str);

			// we ignore the time interval currently
			String sql = "select ST_AsText(start_lonlat) as start_location, ST_AsText(end_lonlat) as end_location from nytaxi "
						+ "where start_time>=" + interval_beg_str 
						+ " and start_time<" + interval_end_str
						+";";
				
				ResultSet rs = statement.executeQuery(sql);
				
				System.out.println(interval_beg_str);

				//Go over each trip record
				while (rs.next()) {
										
					String orig_str = rs.getString(1);
					String dest_str = rs.getString(2);
					
					Point_latlon orig_point_latlon = convertWKT2PointLatlon(orig_str);
					Point_latlon dest_point_latlon = convertWKT2PointLatlon(dest_str);
					
					Point orig_point = convert.convertLatLonToUTM(orig_point_latlon);
					Point dest_point = convert.convertLatLonToUTM(dest_point_latlon);
					
					int orig_grid_index = locatePointInGridNeighbor(orig_point);
					int dest_grid_index = locatePointInGridNeighbor(dest_point);
					//int hour_index = interval_beg.get(Calendar.HOUR_OF_DAY);
					if(orig_grid_index != -1 && dest_grid_index != -1)
					{
						int orig_cid = grid2cluster.get(orig_grid_index);
						int dest_cid = grid2cluster.get(dest_grid_index);
						hour_flow[orig_cid][dest_cid]++;
					}
				}
				
				int hour_index = interval_beg.get(Calendar.HOUR_OF_DAY);
				
				for(int i=0;i<nneigh;i++)//original
				{
					int sum = 0;
					for(int k=0;k<nneigh;k++)
						sum += hour_flow[i][k];
					
					String record =sum + " " + Integer.toString(i+1) + " " + Integer.toString(hour_index+1);
					//ofeature_Writer.write(feature_line + '\n');
					o_flow_Writer.write(record+ '\n');
				}
				
				ctime.add(Calendar.HOUR, 1);
			}
		
		o_flow_Writer.close();
		//ofeature_Writer.close();
		//dest_featureWriter.close();
		//dest_valWriter.close();
		System.out.println("Finish Gaussian Training Data");
	}
	
	
	public boolean generateTheTensor(String filename, Calendar btime, Calendar etime)
			throws SQLException, IOException {
		
		//Step0: initialize the tensor
		mtensor = new int[tensor_time_period][][];
		
		for(int k=0;k<tensor_time_period;k++){//for each time mode
			mtensor[k] = new int[ncluster][];
			for(int i=0;i<ncluster;i++)//for each original mode
				mtensor[k][i] = new int[ncluster];
		}
		
		//Step1: Get the trip records from the database
		CoordinateConversion convert = new CoordinateConversion();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		Calendar ctime = (Calendar)btime.clone();
		Calendar pretime = (Calendar)ctime.clone();
		int workday_number = 1;
		
		while(ctime.before(etime)) {
			if(ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SUNDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SATURDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.FRIDAY){
				ctime.add(Calendar.HOUR_OF_DAY, 1);
				continue;
			}
			

			if(pretime.get(Calendar.DAY_OF_MONTH) != ctime.get(Calendar.DAY_OF_MONTH))
				workday_number++;
			//String str_grid_rect = transfer2PostGISPolygon(llpoint, rupoint);
			Calendar interval_beg = (Calendar)ctime.clone();
			Calendar interval_end = (Calendar)ctime.clone();
			interval_end.add(Calendar.HOUR, 1);
			String interval_beg_str = transfer2TimeStampStr(interval_beg,"EDT");
			String interval_end_str = transfer2TimeStampStr(interval_end,"EDT");
			//System.out.println(interval_beg_str);

			// we ignore the time interval currently
			String sql = "select ST_AsText(start_lonlat) as start_location, ST_AsText(end_lonlat) as end_location from nytaxi "
						+ "where start_time>=" + interval_beg_str 
						+ " and start_time<" + interval_end_str
						+";";
				
				ResultSet rs = statement.executeQuery(sql);
				
				//////////////////////// for debug
				int count_whole = 0; // total number of trips
				int count_within = 0; // the number of trips within query neighborhood
				/////////////////////////

				while (rs.next()) {
										
					String orig_str = rs.getString(1);
					String dest_str = rs.getString(2);
					
					Point_latlon orig_point_latlon = convertWKT2PointLatlon(orig_str);
					Point_latlon dest_point_latlon = convertWKT2PointLatlon(dest_str);
					
					Point orig_point = convert.convertLatLonToUTM(orig_point_latlon);
					Point dest_point = convert.convertLatLonToUTM(dest_point_latlon);
					
					int orig_grid_index = locatePointInGridNeighbor(orig_point);
					int dest_grid_index = locatePointInGridNeighbor(dest_point);
					int hour_index = interval_beg.get(Calendar.HOUR_OF_DAY);
					if(orig_grid_index != -1 && dest_grid_index != -1)
					{
						int orig_cid = grid2cluster.get(orig_grid_index);
						int dest_cid = grid2cluster.get(dest_grid_index);
						mtensor[hour_index][orig_cid][dest_cid] += 1;
						count_within+=1;//////////////////////// for debug
					}
					count_whole++;//////////////////////// for debug
				}
				System.out.println("time:"+interval_beg_str+",count_whole:"+count_whole + ",count_within:"+count_within);


				pretime = (Calendar)ctime.clone();
				ctime.add(Calendar.HOUR, 1);
			}
		

		//Step2: Output the tensor to a file
		System.out.println("workday number:"+workday_number);
		FileWriter fwriter = new FileWriter(filename);
		for(int t=0;t<tensor_time_period;t++){//go over times
			//fwriter.append("t\n");
			for(int i=0;i<ncluster;i++){//go over originals
				fwriter.append(Integer.toString(mtensor[t][i][0] / workday_number));
				for(int j=1;j<ncluster;j++){//go over destinations
					fwriter.append(" " + Integer.toString(mtensor[t][i][j]/ workday_number));
				}
				fwriter.append("\n");
			}
		}
		fwriter.close();
		
		System.out.println("Finished writing the tensor file");

		return true;
	}
	
	public static void main(String[] args) throws SQLException, IOException{
		Grid_Mobility gm = new Grid_Mobility("../data/cluster.txt");
//		Polygon poly = new Polygon();
//		poly.addPoint(gm.boundPoints[0]);
//		poly.addPoint(gm.boundPoints[1]);
//		poly.addPoint(gm.boundPoints[2]);
//		poly.addPoint(gm.boundPoints[3]);
//		System.out.println(gm.toPostGISWKT(poly));
		
		
		
		///////////////////Generate the tensor//////////////////////////////
	
		Calendar btime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		btime.set(2013, 8, 3, 0, 0, 0);//the month begins from zero.
		
		Calendar etime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		etime.set(2013, 9, 22, 0, 0, 0);//the month begins from zero.
		
		gm.generateTheTensor("../data/tensor.txt", btime, etime);
		
		/////////////////////////////////////////////////////////////////
		
		
		
		
		/////////////////Generate the Gaussian Process Training Data/////////////
		//gm.gen_gp_inputs("../data/flow_records.txt");
				
		///////////////////////////////////////////////////////////////////////
		
		
		
		
		////////////////////Flow Distribution/////////////////////////////////
		//for outflow, select neighborhood 4, 6, 10, 14, 20, hour 8, 10, 12, 18, 3		
		//gm.flowDistribution(21, 8);
		//System.out.println();
		//gm.flowDistribution(21, 10);		

		//gm.temporalFlow(18);
		//gm.temporalFlow(19);
		//gm.temporalFlow(21);
		//gm.temporalFlow(22);
		
		/////////////////////////////////////////////////////////////////////
		
		
		
		
		
		
		
		///////////////////////////Output the boundary to the file///////////////////
		//int[] a = new int[24];
		//for(int i=0;i<24;i++)
		//	a[i] = i;
		//gm.outputGridsToFile("../data/cluster_boundaries.csv", a);
		
	}
}
