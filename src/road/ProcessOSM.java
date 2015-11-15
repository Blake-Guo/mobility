package road;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.transform.stream.StreamSource;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.Element;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;

import java.util.HashSet;
import java.util.Set;

public class ProcessOSM {
	
	
	/**
	 * Get rid of those ways and corresponding nodes not for vehicle.
	 * @param ifilePath
	 * @param ofilePath
	 * @throws ParserConfigurationException 
	 * @throws IOException 
	 * @throws SAXException 
	 */
	public void processOSMFile(String ifilePath, String ofilePath) {
		try{
			File osmFile = new File(ifilePath);
			DocumentBuilder dBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder();
		
			Document doc = dBuilder.parse(osmFile);
			doc.normalize();
			Node root = doc.getFirstChild();
			System.out.println(root.getNodeName());
			
			Node node_way = root.getFirstChild();
			//remove unsuitable ways.
//			removeRecursively(doc,Node.ELEMENT_NODE,"way");
			while(node_way != null)
			{
				if(node_way.getNodeName().equals("way")){
					if(node_way.getNodeType() == Node.ELEMENT_NODE && isValidWayforVechile(node_way) == false){
						System.out.println(node_way.getNodeName());
						node_way.getParentNode().removeChild(node_way);
					}
				}
				node_way = node_way.getNextSibling();
			}
			
			//write out the xml
			TransformerFactory factory = TransformerFactory.newInstance();
			Transformer transformer = factory.newTransformer();
			DOMSource source = new DOMSource(doc);         
			StreamResult result = new StreamResult(new File(ofilePath)); 
			transformer.transform(source, result);
			

			System.out.println("Done");
		}
		catch( Exception e){
			e.printStackTrace();
		}
	}
	
	public void removeRecursively(Node node, short nodeType, String name) {
		if (node.getNodeType()==nodeType && (name == null || node.getNodeName().equals(name))) {
			if(isValidWayforVechile(node)==false)
				node.getParentNode().removeChild(node);
		}
		else {
			// check the children recursively
			NodeList list = node.getChildNodes();
			for (int i = 0; i < list.getLength(); i++) {
				removeRecursively(list.item(i), nodeType, name);
			}
		}
	}
	
	
	/**
	 * Check whether the given node is a way valid for vechile.
	 * @param node_way
	 * @return
	 */
	public boolean isValidWayforVechile(Node node_way){
		
		Set<String> forbiddenWays = new HashSet<String>();
		
		forbiddenWays.add("service");
		forbiddenWays.add("pedestrian");
		forbiddenWays.add("track");
		forbiddenWays.add("raceway");
		forbiddenWays.add("footway");
		forbiddenWays.add("bridleway");
		forbiddenWays.add("steps");
		forbiddenWays.add("path");
		forbiddenWays.add("cycleway");
		forbiddenWays.add("proposed");
		forbiddenWays.add("bus_stop");
		forbiddenWays.add("crossing");
		forbiddenWays.add("elevator");
		forbiddenWays.add("emergency_access_point");
		forbiddenWays.add("escape");
		forbiddenWays.add("phone");
		forbiddenWays.add("rest_area");
		forbiddenWays.add("traffic_signals");

		
		NodeList list_wc = node_way.getChildNodes();
		boolean isHighway = false;
		
		for(int j=0;j<list_wc.getLength();j++){
			Node node_way_child = list_wc.item(j);
			
			//=====for test====

			//continue;
			//=================
			
			if(node_way_child.getNodeName().equals("tag") && node_way_child.hasAttributes()){
				NamedNodeMap nodeMap = node_way_child.getAttributes();
				if( nodeMap.item(0).getNodeValue().equals("highway")){
					isHighway = true;
					String nodeval = nodeMap.item(1).getNodeValue();
					if(forbiddenWays.contains(nodeval)){
						return false;
					}
					break;
				}
			}
		}
		
		return isHighway;
	}
	
	
	public void printNote(NodeList nodeList) {

	    for (int count = 0; count < nodeList.getLength(); count++) {

		Node tempNode = nodeList.item(count);

		// make sure it's element node.
		if (tempNode.getNodeType() == Node.ELEMENT_NODE) {

			// get node name and value
			System.out.println("\nNode Name =" + tempNode.getNodeName() + " [OPEN]");

			if (tempNode.hasAttributes()) {

				// get attributes names and values
				NamedNodeMap nodeMap = tempNode.getAttributes();

				for (int i = 0; i < nodeMap.getLength(); i++) {

					Node node = nodeMap.item(i);
					System.out.println("attr name : " + node.getNodeName());
					System.out.println("attr value : " + node.getNodeValue());

				}
			}

			if (tempNode.hasChildNodes()) {

				// loop again if has child nodes
				printNote(tempNode.getChildNodes());

			}

			System.out.println("Node Name =" + tempNode.getNodeName() + " [CLOSE]");
		}
	    }


	    }

	
	public static void main(String[] args){
		ProcessOSM posm = new ProcessOSM();
		posm.processOSMFile("data/newyork_map.osm", "data/new_map.osm");
	}
}
