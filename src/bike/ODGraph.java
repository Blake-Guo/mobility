package bike;

import java.io.*;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.*;

import Geo.Point_latlon;
import road.*;



/*
 * Instruction: for building the ODgraph, just initialize the ODGraph with the constructor ODGraph(String tripFile) or call initOD(tripFile) later, 
 * then everything will be taken care of.
 */

public class ODGraph {

	int[][] od_mat;// the OD matrix
	Map<Integer, Station> stations;// the station collection
	Map<Integer, Integer> id2ind;// convert station's id to a corresponding
									// index
	Map<Integer, Integer> ind2id;// convert station's index to its id.

	public ODGraph() {

	}
	

	/**
	 * By passing the path of trip records, this constructor will construct and build everything well including od_mat, stations, id2ind and ind2id.
	 * @param tripFile
	 * @throws IOException
	 */
	public ODGraph(String tripFile) throws IOException {
		initOD(tripFile);
	}
	
	
	/**
	 * By passing the path of trip records, this constructor will construct and build everything well including od_mat, stations, id2ind and ind2id.
	 * @param tripFile
	 * @param STFlag, indicate whether we want the temporal feature vector or the original spatial feature vector.
	 * @throws IOException
	 */
	public ODGraph(String tripFile, boolean STFlag) throws IOException {
		if(STFlag == false){
			initOD(tripFile);
		}
		else{
			initSTOD(tripFile);
		}
	}

	public void initStations_Jason(String http) {

	}

	/**
	 * Parse one line of trip records and convert it into a station.
	 * 
	 * @param strs
	 * @param ssFlag, mark whether we want to parse the start station or not.
	 * @return
	 */
	public Station parse2Station(String[] strs, boolean ssFlag) {

		Station s = new Station();
		
		if (ssFlag) {//start station
			s.id = Integer.valueOf(strs[3]);
			s.address = strs[4];
			s.point_latlon = new Point_latlon(Double.valueOf(strs[5]),
					Double.valueOf(strs[6]));
		}
		else{//end station
			s.id = Integer.valueOf(strs[7]);
			s.address = strs[8];
			s.point_latlon = new Point_latlon(Double.valueOf(strs[9]),
					Double.valueOf(strs[10]));
		}

		CoordinateConversion convert = new CoordinateConversion();

		s.point_xy = convert.convertLatLonToUTM(s.point_latlon.lat,
				s.point_latlon.lon);

		return s;
	}

	
	/**
	 * Build the station map, id2ind and ind2id from trip data 
	 * @param tripFile, the path of the trip data.
	 * @throws IOException
	 */
	public void initStations_TripData(String tripFile) throws IOException {
		stations = new HashMap<Integer, Station>();
		id2ind = new TreeMap<Integer, Integer>();
		ind2id = new TreeMap<Integer, Integer>();

		//1: Read in the trip file.
		BufferedReader br = new BufferedReader(new FileReader(tripFile));

		String line = br.readLine(); // skip the title

		//2: build the station map
		while ((line = br.readLine()) != null) {
			String[] strs = line.split(",");
			for(int i=0;i<strs.length;i++)
				strs[i] = strs[i].substring(1,strs[i].length()-1);

			// 0:"tripduration",1:"starttime",2:"stoptime",3:"start station id",4:"start station name",5:"start station latitude",6:"start station longitude",
			// 7:"end station id",8:"end station name",9:"end station latitude",10:"end station longitude",11:"bikeid",12:"usertype",13:"birth year",14:"gender"
			int sid = Integer.valueOf(strs[3]);
			int eid = Integer.valueOf(strs[7]);
			
			//start station
			if(stations.containsKey(sid) == false){
				Station s = parse2Station(strs,true);
				stations.put(sid, s);
			}
			
			//end station
			if(stations.containsKey(eid) == false){
				Station s = parse2Station(strs,false);
				stations.put(eid, s);
			}

		}
		
		
		//3: initialize the id2ind and ind2id based on the station map
		int ind = 0;
		for(Map.Entry<Integer, Station> entry : stations.entrySet()){
			id2ind.put(entry.getKey(), ind);
			ind2id.put(ind, entry.getKey());
			ind++;
		}

		br.close();
		
		System.out.println("Finish building the station map, id2ind and ind2id from trip data");
	}

	/**
	 * Initialize the OD matrix
	 * 
	 * @param tripFile
	 *            , the path of the trip data file
	 * @throws IOException
	 */
	public void initOD(String tripFile) throws IOException {

		// 0: Initialize the station map
		initStations_TripData(tripFile);

		int nStation = stations.size();
		
		// 1: Construct the OD matrix
		od_mat = new int[nStation][];
		for (int i = 0; i < nStation; i++)
			od_mat[i] = new int[nStation];

		// 2: Read in the trip file and build the OD matrix.
		BufferedReader br = new BufferedReader(new FileReader(tripFile));

		String line = br.readLine(); // skip the title

		while ((line = br.readLine()) != null) {
			String[] strs = line.split(",");
			for(int i=0;i<strs.length;i++)
				strs[i] = strs[i].substring(1,strs[i].length()-1);
			
			// 0:"tripduration",1:"starttime",2:"stoptime",3:"start station id",4:"start station name",5:"start station latitude",6:"start station longitude",
			// 7:"end station id",8:"end station name",9:"end station latitude",10:"end station longitude",11:"bikeid",12:"usertype",13:"birth year",14:"gender"
			
			int sid = Integer.valueOf(strs[3]);
			int eid = Integer.valueOf(strs[7]);
			
			od_mat[id2ind.get(sid)][id2ind.get(eid)]++;
		}

		br.close();

		System.out.println("Finish initializing the OD matrix");
		
		//System.out.println("the number of all trips:" + sumofODTrip());

	}
	
	
	
	
	/**
	 * Initialize the OD matrix with temporal feature, [  [destinations in the morning], [destinations in the afternoon]， [destinations in the evening]  ] 
	 * 
	 * @param tripFile
	 *            , the path of the trip data file
	 * @throws IOException
	 * @throws  
	 */
	public void initSTOD(String tripFile) throws IOException {

		// 0: Initialize the station map
		initStations_TripData(tripFile);

		int nStation = stations.size();
		
		// 1: Construct the temporal OD matrix
		od_mat = new int[nStation][];
		for (int i = 0; i < nStation; i++)
			od_mat[i] = new int[nStation * 3]; // multiple 3 because there are morning, afternoon and evening time periods.

		// 2: Read in the trip file and build the OD matrix.
		BufferedReader br = new BufferedReader(new FileReader(tripFile));

		String line = br.readLine(); // skip the title

		while ((line = br.readLine()) != null) {
			String[] strs = line.split(",");
			for(int i=0;i<strs.length;i++)
				strs[i] = strs[i].substring(1,strs[i].length()-1);//because each string has the symbol of '"', '"'.
			
			// 0:"tripduration",1:"starttime",2:"stoptime",3:"start station id",4:"start station name",5:"start station latitude",6:"start station longitude",
			// 7:"end station id",8:"end station name",9:"end station latitude",10:"end station longitude",11:"bikeid",12:"usertype",13:"birth year",14:"gender"
			
			int sid = Integer.valueOf(strs[3]);
			int eid = Integer.valueOf(strs[7]);
			
			//get the start time, its format is "2014-07-01 00:00:04"
			DateFormat data_format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
			java.util.Date stime = null;
			try {
				stime = data_format.parse(strs[2]);
			} catch (ParseException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			Calendar cal = Calendar.getInstance();
			cal.setTime(stime);
			
			
			int hour = cal.get(Calendar.HOUR_OF_DAY);
			if(hour > 6 && hour < 12){
				od_mat[id2ind.get(sid)][id2ind.get(eid)]++;
			}
			else if(hour >= 12 && hour < 17){
				od_mat[id2ind.get(sid)][id2ind.get(eid) + nStation]++;
			}
			else if(hour >= 17 && hour <= 23){
				od_mat[id2ind.get(sid)][id2ind.get(eid) + 2*nStation]++;
			}
			
		}

		br.close();

		System.out.println("Finish initializing the Spatio-Temporal OD matrix");
		
	}
	
	
	
	public int sumofODTrip(){
		int sum=0;
		for(int i=0;i<od_mat.length;i++)
			for(int j=0;j<od_mat[i].length;j++)
				sum += od_mat[i][j];
		
		return sum;
	}
	
	
	
	/**
	 * Return the top k distance within station sId(index, not ID).
	 * @param sInd
	 * @param k
	 * @return
	 */
	public PriorityQueue<Double> topKLongestLen(int sInd, int k){
		
		PriorityQueue<Double> pq = new PriorityQueue<Double>(k, Collections.reverseOrder());
		int id = ind2id.get(sInd);
		Station s = stations.get(id);
		for(Map.Entry<Integer, Station> entry : stations.entrySet()){
			if(entry.getKey() == id)
				continue;
			
			double dis = s.distance_man(entry.getValue());
			
			if(pq.size() == k)
			{
				if(pq.peek() > dis){
					pq.poll();
					pq.add(dis);
				}
			}
			else{
				pq.add(dis);
			}
		}
		
		return pq;
	}
	
	
	/**
	 * Output the station map to a file with title "station_ind", "latitude", "longitude"
	 * @param filepath
	 * @throws IOException
	 */
	public void outputStationMap(String filepath)throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(filepath));
		String line = "station_ind" + "\t" + "latitude" + "\t" + "longitude";
		bw.write(line);
		bw.newLine();
		
		
		for(int i=0;i<stations.size();i++){
			int id = ind2id.get(i);//from index to ID
			
			Station s = stations.get(id);
			
			line = Integer.toString(i) + "\t" + Double.toString(s.point_latlon.lat) + "\t" + Double.toString(s.point_latlon.lon);
			bw.write(line);
			bw.newLine();
			bw.flush();
		}
		
		bw.close();
	}
	
	
	public void outputTripsRepWithInd(String inTripPath, String outTripPath) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(inTripPath));
		BufferedWriter bw = new BufferedWriter(new FileWriter(outTripPath));
		
		String line;
		
		//title;
		line = br.readLine();
		String[] strs = line.split(",");
		strs[3] = "start station ind";
		strs[7] = "end station ind";
		line = concatStrings(strs,',');
		bw.write(line);
		bw.newLine();
		
		while((line = br.readLine()) != null){
			strs = line.split(",");
			int s_ind = id2ind.get(Integer.valueOf(strs[3].substring(1,strs[3].length()-1)));
			int e_ind = id2ind.get(Integer.valueOf(strs[7].substring(1,strs[7].length()-1)));
			
			strs[3] = '"' + Integer.toString(s_ind) + '"';
			strs[7] = '"' + Integer.toString(e_ind)+ '"';
			
			line = concatStrings(strs,',');
			bw.write(line);
			bw.newLine();
			bw.flush();
		}
		
		bw.close();
		br.close();
		
		
		System.out.println("Finish output the trips replaced with index");
		return;
	}
	
	
	public String concatStrings(String[] strs, char sep){
		String line = strs[0];
		
		for(int i=1;i<strs.length;i++)
			line += (sep + strs[i]);
		
		return line;
	}
	

	public static void main(String[] args) throws IOException {
		ODGraph b_graph = new ODGraph("data/2014-07-trip data.csv");
		//b_graph.outputTripsRepWithInd("data/2014-07-trip data.csv", "bike_trips.csv");
		
		System.out.println(b_graph.topKLongestLen(154, 20));
		//b_graph.outputStationMap("bike_stations.csv");
		
	}
}
