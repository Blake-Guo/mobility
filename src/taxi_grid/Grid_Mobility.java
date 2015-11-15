package taxi_grid;

import utility.ConnectDB;

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
import java.util.Locale;
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
	double grid_width;// this is the width of a grid, and the unit metric is kilometer
	double grid_height;// this is the height of a grid, and the unit metric kilometer
	double boundary_width;// the width of the boundary
	double boundary_height;// the height of the boundary
	Polygon[] grid_boundaries;
	int grid_xnumber;
	int grid_ynumber;
	int grid_number;
	int[][][] mtensor;//the mode of mtensor is reverse of the real mobility tensor. Here mode 1: time, mode 2: original, mode 3: destination
	final int tensor_time_period = 24;
	
	public Grid_Mobility(){
		//The boundary of Manhatten, (-74.029999,40.703546), (-73.974380,40.703546), (-73.924942,40.791459),(-73.98056,40.791459)
		CoordinateConversion convert = new CoordinateConversion();
		boundPoints = new Point[4];
		boundPoints[0] = convert.convertLatLonToUTM(40.703546, -74.029999);
		boundPoints[1] = convert.convertLatLonToUTM(40.703546, -73.974380);
		boundPoints[2] = convert.convertLatLonToUTM(40.791459, -73.924942);
		boundPoints[3] = convert.convertLatLonToUTM(40.791459, -73.98056);
		
		//1: initialize the grids
		grid_width = 800;//1kilometer
		grid_height = 800;//1 kilometer
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
		
		//2: initialize the tensor
		mtensor = new int[tensor_time_period][][];
		
		for(int k=0;k<tensor_time_period;k++){//for each time mode
			mtensor[k] = new int[grid_number][];
			for(int i=0;i<grid_number;i++)//for each original mode
				mtensor[k][i] = new int[grid_number];
		}
		
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
//		int i = index / grid_ynumber;
//		int j = index % grid_ynumber;
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
	
	
	public boolean generateTheTensor(String filename, Calendar btime, Calendar etime)
			throws SQLException, IOException {
		
		//Step1: Get the trip records from the database
		CoordinateConversion convert = new CoordinateConversion();
		Connection conn = ConnectDB.ConnectDatabase("taxidb");
		Statement statement = conn.createStatement();
		Calendar ctime = (Calendar)btime.clone();
		while(ctime.before(etime)) {
			if(ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SUNDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.SATURDAY || ctime.get(Calendar.DAY_OF_WEEK) == Calendar.FRIDAY){
				ctime.add(Calendar.HOUR_OF_DAY, 1);
				continue;
			}

			//String str_grid_rect = transfer2PostGISPolygon(llpoint, rupoint);
			Calendar interval_beg = (Calendar)ctime.clone();
			Calendar interval_end = (Calendar)ctime.clone();
			interval_end.add(Calendar.HOUR, 1);
			String interval_beg_str = transfer2TimeStampStr(interval_beg,"EDT");
			String interval_end_str = transfer2TimeStampStr(interval_end,"EDT");
			System.out.println(interval_beg_str);

			// we ignore the time interval currently
			String sql = "select ST_AsText(start_lonlat) as start_location, ST_AsText(end_lonlat) as end_location from nytaxi "
						+ "where start_time>=" + interval_beg_str 
						+ " and start_time<" + interval_end_str
						+";";
				
				ResultSet rs = statement.executeQuery(sql);

				while (rs.next()) {
					
					String orig_str = rs.getString(1);
					String dest_str = rs.getString(2);
					
					Point_latlon orig_point_latlon = convertWKT2PointLatlon(orig_str);
					Point_latlon dest_point_latlon = convertWKT2PointLatlon(dest_str);
					
					Point orig_point = convert.convertLatLonToUTM(orig_point_latlon);
					Point dest_point = convert.convertLatLonToUTM(dest_point_latlon);
					
					int orig_grid_index = locatePointInGridNeighbor(orig_point);
					int dest_grid_index = locatePointInGridNeighbor(dest_point);
					int hour_index = interval_beg.get(Calendar.HOUR);
					if(orig_grid_index != -1 && dest_grid_index != -1)
						mtensor[hour_index][orig_grid_index][dest_grid_index] += 1;
				}

				ctime.add(Calendar.HOUR, 1);
			}

		//Step2: Output the tensor to a file
		FileWriter fwriter = new FileWriter(filename);
		for(int t=0;t<tensor_time_period;t++){//go over times
			fwriter.append("t\n");
			for(int i=0;i<grid_number;i++){//go over originals
				fwriter.append(Integer.toString(mtensor[t][i][0]) + " ");
				for(int j=1;j<grid_number;j++){//go over destinations
					fwriter.append(Integer.toString(mtensor[t][i][j]) + " ");
				}
				fwriter.append("\n");
			}
		}
		fwriter.close();
		
		System.out.println("Finished writing the tensor file");

		return true;
	}
	
	public static void main(String[] args) throws SQLException, IOException{
		Grid_Mobility gm = new Grid_Mobility();
		
		Calendar btime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		btime.set(2013, 8, 3, 0, 0, 0);//the month begins from zero.
		//System.out.println(btime.getTimeZone().toString());
		
		Calendar etime = new GregorianCalendar(TimeZone.getTimeZone("GMT-4"));
		etime.set(2013, 9, 31, 0, 0, 0);//the month begins from zero.
		
		//SimpleDateFormat f = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		//f.setTimeZone(TimeZone.getTimeZone("GMT-4"));
		//String timestamp1 = f.format(btime.getTime());
		//System.out.println(timestamp1);
		
		
		gm.generateTheTensor("data/tensor.txt", btime, etime);
	}
}
