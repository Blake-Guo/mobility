package Geo;

public class IntersectTest {
	
	/**
	 * Check whether segment(p1,p2) intersect with segment(p3,p4) or not
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @return
	 */
	public static boolean segmentInterset(Point p1, Point p2, Point p3, Point p4){
		double d1 = direction(p3,p4,p1);
		double d2 = direction(p3,p4,p2);
		double d3 = direction(p1,p2,p3);
		double d4 = direction(p1,p2,p4);
		
		if( (d1*d2<0) && (d3*d4<0) )
			return true;
		else if(d1 == 0 && onSegmentRange(p3,p4,p1))
			return true;
		else if(d2 == 0 && onSegmentRange(p3,p4,p2))
			return true;
		else if(d3 == 0 && onSegmentRange(p1,p2,p3))
			return true;
		else if(d4 == 0 && onSegmentRange(p1,p2,p4))
			return true;
		return false;
	}
	
	/**
	 * Use equation to see if ray (p with direction dir) intersects segment (a,b)
	 * This equation comes from https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
	 * @param o
	 * @param dir
	 * @param a
	 * @param b
	 * @return
	 */
	public static boolean rayIntersectSegment(Point o, Point dir, Point a, Point b){
		Point v1 = new Point(o.x - a.x, o.y - a.y);
		Point v2 = new Point(b.x - a.x, b.y - a.y);
		Point v3 = new Point(-dir.y, dir.x);
		
		double t1 = crossProduct(v2,v1) / dotProduct(v2,v3);
		double t2 = dotProduct(v1,v3) / dotProduct(v2,v3);
		
		if(t1 >=0 && t2>=0 && t2<=1)
			return true;
		return false;
	}
	
	
	/**
	 * Just the cross product
	 * @param v1
	 * @param v2
	 * @return
	 */
	public static double crossProduct(Point v1, Point v2){
		return v1.x * v2.y - v1.y * v2.x;
	}
	
	public static double dotProduct(Point v1, Point v2){
		return v1.x * v2.x + v1.y * v2.y;
	}
	
	
	private static double direction(Point p1, Point p2, Point p3){
		Point v1 = new Point(p3.x - p1.x, p3.y - p1.y);
		Point v2 = new Point(p2.x - p1.x, p2.y - p1.y);
		return crossProduct(v1,v2);
	}
	
	/**
	 * Check whether point p is within the range of segment (p1,p2)
	 * @param p1
	 * @param p2
	 * @param p
	 * @return
	 */
	private static boolean onSegmentRange(Point p1, Point p2, Point p){
		boolean flag1 = (Math.min(p1.x, p2.x) <= p.x) && (p.x <= Math.max(p1.x, p2.x));
		boolean flag2 = (Math.min(p1.y, p2.y) <= p.y) && (p.y <= Math.max(p1.y, p2.y));
		return flag1 && flag2;
	}
	
	
	
}
