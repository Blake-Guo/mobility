package road;
import java.util.*;


public class Edge {
	public int edgeId;
	public String snodeId;
	public String enodeId;
	
	public List<Point> points;//the points to represent the shape of the edge including the start node and end node.
	//public List<Integer> adjEdgeIds;//the id of the edges attached to this edge.
	
	public Edge(){
		edgeId = -1;
		snodeId = "";
		enodeId = "";
		
		points = new ArrayList<Point>();
		//adjEdgeIds = new LinkedList<Integer>();
	}

}
