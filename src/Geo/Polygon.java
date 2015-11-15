package Geo;
import java.util.ArrayList;


public class Polygon {
	public ArrayList<Point> points;//currently the last point is not the same as first point.
	
	public Polygon(){
		points = new ArrayList<Point>();
	}
	
	public Polygon(ArrayList<Point> list){
		points = new ArrayList<Point>(list);
	}
	
	public void addPoint(Point p){
		points.add(p);
	}
	
	public boolean contain(Point p){
		boolean flag1 = contains1(p);
		boolean flag2 = contains2(p);
		
		if(flag1 != flag2){
			try {
				throw new Exception("The point within polygon results are different");
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				//print the polygon
				System.out.println("The polygon is:");
				for(int i=0;i<points.size();i++){
					System.out.println(points.get(i).x + "," + points.get(i).y);
				}
				//print the point
				System.out.println("The point to be tested is:" + p.x + "," + p.y);
			}
		}
		return flag1;
	}
	
	/**
	 * Use the two segment intersection to test whether the given point p is within the polygon
	 * @param p
	 * @return
	 */
	private boolean contains1(Point p){
		Point inf_point = new Point(p.x, 18840000);
		int count = 0;
		for(int i=1;i<points.size();i++){
			if(IntersectTest.segmentInterset(points.get(i-1), points.get(i), p, inf_point))
				count++;
		}
		int n = points.size();
		if(IntersectTest.segmentInterset(points.get(n-1), points.get(0), p, inf_point))
			count++;
		if(count % 2 == 0)
			return false;
		return true;
	}
	
	/**
	 * Use the ray intersect with segment equation.
	 * @param p
	 * @return
	 */
	private boolean contains2(Point p){
		Point dir = new Point(0.5, 0.8660254);
		int count = 0;
		for(int i=1;i<points.size();i++){
			if(IntersectTest.rayIntersectSegment(p,dir, points.get(i-1), points.get(i)) )
				count++;
		}
		int n = points.size();
		if(IntersectTest.rayIntersectSegment(p,dir, points.get(n-1), points.get(0)))
			count++;
		if(count % 2 == 0)
			return false;
		return true;
	}
}
