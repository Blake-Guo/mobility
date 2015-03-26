package road;

import java.util.*;
import java.awt.geom.*;

//^ y height
//|
//|
//|
//---------->x width
public class GridIndex {
	public Point swpoint;
	public Point nepoint;
	public int nrow = 50;//row starts from 0
	public int ncol = 50;//col starts from 0.
	public double gridw;
	public double gridh;

	public LinkedList<Integer>[][] grids; //grids[i][j] store the id of edges which intercept with grids[i][j]

	
	public LinkedList<Integer> getNearbyEdgeIds(Point p){
	    int rowId, colId;
	    rowId =(int) ((p.y - swpoint.y) / gridh);
	    colId =(int) ((p.x - swpoint.x) / gridw);
	    rowId = (rowId < 0) ? 0 : rowId;
	    rowId = (rowId >= nrow) ? (nrow-1) : rowId;
	    
	    colId = (colId < 0) ? 0 : colId;
	    colId = (colId >= ncol) ? (ncol-1) : colId;
	    
	    return grids[rowId][colId];
	}
	
	/**
	 * Initialize the grid indexing based on the given map
	 * @param map
	 */
	public void initGridTree(RoadNetwork map) {
		
		swpoint = map.point_bl;
		nepoint = map.point_tr;
		
		gridw = (map.point_tr.x - map.point_bl.x) / ncol;
		gridh = (map.point_tr.y - map.point_bl.y) / nrow;
		
		// initialize the grids based on the given map.
		grids = (LinkedList<Integer>[][]) new LinkedList<?>[nrow][];
		
		for(int i=0;i<nrow;i++)
		{
			grids[i] = (LinkedList<Integer>[]) new LinkedList<?>[ncol];
			for(int j =0;j<ncol;j++)
				grids[i][j] =  new LinkedList<Integer>();
		}
		
		//scan each edge
		for(Map.Entry<Integer, Edge> entry : map.edges.entrySet()){
			Edge e = entry.getValue();
			Set<Integer> vGrids = new HashSet<Integer>();
			
			
			//scan each point of the edge
			List<Point> points = e.points;
			for(int i=1;i<points.size();i++){
				int brow = (int) ( (points.get(i-1).y - map.point_bl.y) / gridh);
				int trow = (int) ( (points.get(i).y - map.point_bl.y) / gridh);
				if(brow == nrow)
					brow--;
				if(trow == nrow)
					trow--;
				if(brow > trow)
				{
					int tmp = brow;
					brow = trow;
					trow = tmp;
				}
				
				int lcol = (int) ( (points.get(i-1).x - map.point_bl.x ) / gridw);
				int rcol = (int) ( (points.get(i).x - map.point_bl.x) / gridw);
				if(lcol == ncol)
					lcol--;
				if(rcol == ncol)
					rcol--;
				if(lcol > rcol)
				{
					int tmp = lcol;
					lcol = rcol;
					rcol = tmp;
				}
				
				
				for(int r=brow;r<=trow;r++)
				for(int c = lcol;c<=rcol;c++)
				{
					int hashval = r * ncol + c;
					if(vGrids.contains(hashval))
						continue;
					
					if(gridIntersect(r,c,points.get(i-1),points.get(i))){
						vGrids.add(hashval);
						grids[r][c].add(e.edgeId);
					}
				}
				
			}
			
			vGrids.clear();
		}
	}
	
	
	/**
	 * Test whether the given grid intersect with the given segment.
	 * @param row
	 * @param col
	 * @param p1
	 * @param p2
	 * @return
	 */
	public boolean gridIntersect(int row, int col, Point p1, Point p2){
		double lx = col * gridw + swpoint.x;
		double by = row * gridh + swpoint.y;
		
		 Rectangle2D grid = new Rectangle2D.Double(lx,by,gridw,gridh);//x,y,w,h
		 return grid.intersectsLine(p1.x, p1.y, p2.x, p2.y);
	}
	
	
}
