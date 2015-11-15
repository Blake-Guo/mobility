package road;

import java.util.*;
import java.awt.geom.Line2D;
import java.io.*;

import Geo.Point;
import Geo.Point_latlon;


public class RoadNetwork {
	public Map<String, Node> nodes;
	public Map<Integer, Edge> edges;
	public GridIndex gindex;
	Point point_bl;
	Point point_tr;
	
	public RoadNetwork(){
		nodes = new HashMap<String, Node>();
		edges = new HashMap<Integer, Edge>();
		gindex = new GridIndex();
		point_bl = new Point(Double.MAX_VALUE,Double.MAX_VALUE);
		point_tr = new Point(Long.MIN_VALUE, Long.MIN_VALUE);
	}
	
	
	public RoadNetwork(String nodeFile, String edgeFile) throws IOException{
		
		nodes = new HashMap<String, Node>();
		edges = new HashMap<Integer, Edge>();
		gindex = new GridIndex();
		point_bl = new Point(Double.MAX_VALUE,Double.MAX_VALUE);
		point_tr = new Point(Long.MIN_VALUE, Long.MIN_VALUE);
		
		initNodes(nodeFile);
		
		initEdges(edgeFile);
	}
	
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
			nodes.get(e.enodeId).adjEdgeIds.add(e.edgeId);//3390355186
			
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
		
		//Update the adjacent edge for all the edges
		for(Map.Entry<Integer, Edge> entry : edges.entrySet()){
			Edge e = entry.getValue();
			
			e.adjEdgeIds.addAll( nodes.get(e.snodeId).adjEdgeIds);
			e.adjEdgeIds.addAll( nodes.get(e.enodeId).adjEdgeIds);
		}
		
		System.out.println("There are " + edges.size() + " edges");

	}
	
	
	/**
	 * Get the neighbor of edge eid within k levels.
	 * @param k
	 * @param eid
	 * @return
	 */
	public Set<Integer> getEdgeKNeighbor(int k, int eid){
		
		Queue<Integer> queue = new LinkedList<Integer>();
		Set<Integer> neighbors = new HashSet<Integer>();
		
		Node snode = nodes.get(edges.get(eid).snodeId);
		Node enode = nodes.get(edges.get(eid).enodeId);

		
		queue.add(eid);
	
		int l = 0;
		
		while(l<=k){
			
			Queue<Integer> newQueue = new LinkedList<Integer>();
			
			while(queue.isEmpty()==false){
				
				int neigh_eid = queue.poll();
				if(neighbors.contains(neigh_eid))
					continue;
				neighbors.add(neigh_eid);

				
				snode = nodes.get(edges.get(neigh_eid).snodeId);
				enode = nodes.get(edges.get(neigh_eid).enodeId);
				
				for(Integer tid : snode.adjEdgeIds)
					if(neighbors.contains(tid) == false)
						newQueue.add(tid);
				for(Integer tid : enode.adjEdgeIds)
					if(neighbors.contains(tid) == false)
						newQueue.add(tid);
			}
			
			queue = newQueue;
			
			l++;
		}
		
		return neighbors;
	}
	
	/**
	 * Check whether edge eid0 is within the k neighborhood of edge eid0.
	 * @param k
	 * @param eid0
	 * @param eid1
	 * @return
	 */
	public boolean isEdgeKNeighbor(int k, int eid0, int eid1){
		Queue<Integer> queue = new LinkedList<Integer>();
		Set<Integer> neighbors = new HashSet<Integer>();
		
		Node snode = nodes.get(edges.get(eid0).snodeId);
		Node enode = nodes.get(edges.get(eid0).enodeId);

		
		queue.add(eid0);
		int l = 0;
		
		while(l<=k){
			
			Queue<Integer> newQueue = new LinkedList<Integer>();
			
			while(queue.isEmpty()==false){
				
				int neigh_eid = queue.poll();
				if(neighbors.contains(neigh_eid))
					continue;
				neighbors.add(neigh_eid);
				
				if(neigh_eid == eid1)
					return true;

				
				snode = nodes.get(edges.get(neigh_eid).snodeId);
				enode = nodes.get(edges.get(neigh_eid).enodeId);
				
				for(Integer tid : snode.adjEdgeIds)
					if(neighbors.contains(tid) == false)
						newQueue.add(tid);
					
				for(Integer tid : enode.adjEdgeIds)
					if(neighbors.contains(tid) == false)
						newQueue.add(tid);
				
			}
			
			queue = newQueue;

			l++;
		}
		
		return false;
	}
	
	/**
	 * Clip the given map and generate the new node and edge file.
	 * @param nodeInFile
	 * @param edgeInFile
	 * @param nodeOutFile
	 * @param edgeOutFile
	 * @throws IOException 
	 */
	public void clipMap(String nodeInFile, String edgeInFile, String nodeOutFile, String edgeOutFile) throws IOException{
		CoordinateConversion cconvert = new CoordinateConversion();

		//The clip rectangle, the left line and the right line
		Point l_point0 = cconvert.convertLatLonToUTM(40.693422, -74.043882);
		Point l_point1 = cconvert.convertLatLonToUTM(40.880062, -73.938139);
		
		Point r_point0 = cconvert.convertLatLonToUTM(40.576435, -73.708544);
		Point r_point1 = cconvert.convertLatLonToUTM(40.856891, -73.708544);
		
		
		BufferedReader br = new BufferedReader(new FileReader(edgeInFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(edgeOutFile));
		Set<String> withinNodeIds = new HashSet<String>();
		//"edge_id","source","target","length","car","car reverse","bike","bike reverse","foot","WKT"
		String line = "";
		line = br.readLine();//scan the titles
		bw.write(line + '\n');
		int edgeId = 0;
		
		while( (line = br.readLine()) != null){
			
			String[] strs = line.split(",");
			
			//store the line which skips the edge_id
			String tmpLine = "";
			for(int i=1;i<strs.length;i++){
				tmpLine+= ("," + strs[i]);
			}
			
			//0: edge_id
			
			//1: source node id
			String snodeId = strs[1];
			
			//2: target node id
			String enodeId = strs[2];
			
			//3: length, 4:car, 5:car reverse, 6: bike, 7: bike reverse, 8: foot
			
			int npoint = strs.length - 9;
			
			//9: LINESTRING(-74.4388088 40.8320984
			strs[9] = strs[9].substring(11);
			
			//the last one: -74.4388088 40.8320984)
			strs[strs.length-1] = strs[strs.length-1].substring(0, strs[strs.length-1].length()-1);
			
			
			//scan all the points
			int outsideClip = 0;
			for(int i=0;i<npoint;i++){
				String[] lonlat = strs[9+i].split(" ");
				Point_latlon p_lonlat = new Point_latlon(Double.valueOf(lonlat[1]),Double.valueOf(lonlat[0]));//longitude, latitude
				//convert the latitude and longitude to the UTM
				Point p = cconvert.convertLatLonToUTM(p_lonlat.lat, p_lonlat.lon);
				
				if( point2SideOfLine(l_point0,l_point1,p) > 0 ||  point2SideOfLine(r_point0,r_point1,p) < 0){
					outsideClip++;
				}
				
			}
			
			if(outsideClip < npoint)
			{
				tmpLine = edgeId + tmpLine;
				bw.write(tmpLine + '\n');
				edgeId++;
				
				withinNodeIds.add(snodeId);
				withinNodeIds.add(enodeId);
			}
			
			//finish process one edge.
		}
		
		br.close();
		bw.close();
		
		fileterNodes(nodeInFile, nodeOutFile, withinNodeIds);
		
		System.out.println("Finish clipping the map");
	}
	
	/**
	 * Judge which side of the line the point locates, return (p1.x-p0.x) * (p.y-p0.y) - (p.x-p0.x) * (p1.y-p0.y)
	 * @param line_p0
	 * @param line_p1
	 * @param point
	 * @return
	 */
	public double point2SideOfLine(Point line_p0, Point line_p1, Point point){
		double x1 = line_p1.x - line_p0.x;
		double y1 = line_p1.y - line_p0.y;
		
		double x2 = point.x - line_p0.x;
		double y2 = point.y - line_p0.y;
		
		return (x1 * y2) - (x2 * y1);
	}
	
	/**
	 * filter those nodes that has no attached edges.
	 * @param nodeInFile
	 * @param nodeOutFile
	 * @param withinNodeIds
	 */
	public void fileterNodes(String nodeInFile, String nodeOutFile, Set<String> withinNodeIds)throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(nodeOutFile));
		BufferedReader br = new BufferedReader(new FileReader(nodeInFile));

		String line = "";
		line = br.readLine();//scan the tiles:"node_id","longitude","latitude"
		bw.write(line + '\n');
		
		while( (line = br.readLine()) != null){
			String[] strs = line.split(",");
			if(withinNodeIds.contains(strs[0])){
				bw.write(line+ '\n');
			}
		}
		
		bw.close();
		br.close();
		
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
		rnetwork.initMap("data/cnodes.csv", "data/cedges.csv");
		Set<Integer> neighbors = rnetwork.getEdgeKNeighbor(1, 157620);
		for(Integer id : neighbors){
			System.out.print(id + " ");
		}
		//rnetwork.clipMap("data/nodes.csv", "data/edges.csv", "data/cnodes.csv", "data/cedges.csv");
	}

}
