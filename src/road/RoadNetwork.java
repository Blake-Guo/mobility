package road;

import java.util.*;
import java.awt.geom.Line2D;
import java.io.*;


public class RoadNetwork {
	public Map<String, Node> nodes = new HashMap<String, Node>();
	public Map<Integer, Edge> edges = new HashMap<Integer, Edge>();
	public GridIndex gindex = new GridIndex();
	Point point_bl = new Point(Double.MAX_VALUE,Double.MAX_VALUE);
	Point point_tr = new Point(Long.MIN_VALUE, Long.MIN_VALUE);
	
	/**
	 * 
	 * @param nodeFilePath
	 */
	public void initNodes(String nodeFilePath)throws IOException{
		CoordinateConversion cconvert = new CoordinateConversion();
		BufferedReader br = new BufferedReader(new FileReader(nodeFilePath));
		br.readLine();//scan the tiles:"node_id","longitude","latitude"
		
		String line = "";//3390355186,-73.8755883,40.8591258
		while( (line = br.readLine()) != null){
			Node node = new Node();
			String[] strs = line.split(",");
			//0: node id
			//node.nodeId = strs[0];
			
			Point_latlon point_latlon = new Point_latlon();
			//1: longitude
			point_latlon.lon = Double.valueOf(strs[1]);
			//2: latitude
			point_latlon.lat = Double.valueOf(strs[2]);
			
			//convert the latitude and longitude to the UTM
			Point point_xy = cconvert.convertLatLonToUTM(point_latlon.lat, point_latlon.lon);
			node.x = point_xy.x;
			node.y = point_xy.y;
			
			//Update the bounding box
			if(node.x < point_bl.x)
				point_bl.x = node.x;
			if(node.x > point_tr.x)
				point_tr.x = node.x;
			if(node.y < point_bl.y)
				point_bl.y = node.y;
			if(node.y > point_tr.y)
				point_tr.y = node.y;
			
			nodes.put(strs[0], node);
		}
		br.close();
		System.out.println("There are " + nodes.size() + " nodes");
	}
	
	/**
	 * Read the edge file and initialize the map's edges based on that. 
	 * @param edgeFilePath
	 * @throws FileNotFoundException 
	 */
	public void initEdges(String edgeFilePath) throws IOException{
		
		CoordinateConversion cconvert = new CoordinateConversion();
		
		BufferedReader br = new BufferedReader(new FileReader(edgeFilePath));
		//"edge_id","source","target","length","car","car reverse","bike","bike reverse","foot","WKT"
		br.readLine();//scan the titles
		
		String line = "";
		while( (line = br.readLine()) != null){
			
			Edge e = new Edge();
			
			String[] strs = line.split(",");
			
			//0: edge_id
			e.edgeId = Integer.valueOf(strs[0]);
			
			//1: source node id
			e.snodeId = strs[1];
			nodes.get(e.snodeId).adjEdgeIds.add(e.edgeId);
			
			//2: target node id
			e.enodeId = strs[2];
			nodes.get(e.enodeId).adjEdgeIds.add(e.edgeId);
			
			//3: length, 4:car, 5:car reverse, 6: bike, 7: bike reverse, 8: foot
			
			int npoint = strs.length - 9;
			
			//9: LINESTRING(-74.4388088 40.8320984
			strs[9] = strs[9].substring(11);
			
			//the last one: -74.4388088 40.8320984)
			strs[strs.length-1] = strs[strs.length-1].substring(0, strs[strs.length-1].length()-1);
			
			
			//scan all the points
			for(int i=0;i<npoint;i++){
				String[] lonlat = strs[9+i].split(" ");
				Point_latlon p_lonlat = new Point_latlon(Double.valueOf(lonlat[1]),Double.valueOf(lonlat[0]));//longitude, latitude
				//convert the latitude and longitude to the UTM
				Point p = cconvert.convertLatLonToUTM(p_lonlat.lat, p_lonlat.lon);
				
				e.points.add(p);
				
				
				//Update the bounding box
				if(p.x < point_bl.x)
					point_bl.x = p.x;
				if(p.x > point_tr.x)
					point_tr.x = p.x;
				if(p.y < point_bl.y)
					point_bl.y = p.y;
				if(p.y > point_tr.y)
					point_tr.y = p.y;
			}
			
			
			//finish process one edge.
			edges.put(e.edgeId, e);
		}
		
		br.close();
		
		System.out.println("There are " + edges.size() + " edges");

	}
	
	/**
	 * Read in the node file and edge file, then initialize the map.
	 * @param nodeFilePath
	 * @param edgeFilePath
	 */
	public void initMap(String nodeFilePath, String edgeFilePath)throws IOException{
		
		//Initialize the nodes and edges
		initNodes(nodeFilePath);
		
		initEdges(edgeFilePath);
		
		System.out.println("Finish generating the road network");
				
		//Generate the grid indexing.
		gindex.initGridTree(this);
		
		System.out.println("Finish generating the grid index");
	}
	
	
	
	/**
	 * Compute the distance from point p to the segment (lp1,lp2)
	 * @param p
	 * @param lp1
	 * @param lp2
	 * @return
	 */
	public double point2SegDist(Point p, Point lp1, Point lp2){
		Line2D line = new Line2D.Double(lp1.x, lp1.y, lp2.x, lp2.y);
		return line.ptSegDist(p.x, p.y);
	}
	
	
	/**
	 * Given one point, use the grid-index to locate which grid it belongs to,
	 * then search and return the id of edge which the point is closest to.
	 * 
	 * @param point
	 * @return
	 */
	public int projectPointToEdge(Point point) {
		
		List<Integer> edgeids = gindex.getNearbyEdgeIds(point);
		double minDis = Double.MAX_VALUE;
		int minEid = -1;
		
		for(Integer eid : edgeids){
			
			List<Point> points = edges.get(eid).points;
			
			for(int i=1;i<points.size();i++){
				double tdis = point2SegDist(point, points.get(i-1), points.get(i));
				if(tdis < minDis){
					minDis = tdis;
					minEid = eid;
				}
			}
		}
		
		return minEid;
	}
	
	
	public static void main(String[] args)throws IOException{
		RoadNetwork rnetwork = new RoadNetwork();
		rnetwork.initMap("data/nodes.csv", "data/edges.csv");
	}

}
