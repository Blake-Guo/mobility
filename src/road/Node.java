package road;
import java.util.*;

public class Node {
	//public String nodeId;
	public double x;
	public double y;
	public List<Integer> adjEdgeIds; // the edges attached to this node.
	
	public Node(){
		//nodeId = "";
		x = 0;
		y = 0;
		adjEdgeIds = new LinkedList<Integer>();
	}
}
