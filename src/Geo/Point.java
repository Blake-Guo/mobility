package Geo;

public class Point {
	public double x;
	public double y;
	
	public Point(){
		x = 0;
		y = 0;
	}
	
	public Point(double x, double y){
		this.x = x;
		this.y = y;
	}
	
	public double distance_euc(Point p){
		double dis = Math.pow(p.x-x, 2) + Math.pow(p.y-y,2);
		dis = Math.sqrt(dis);
		return dis;
	}
	
	public double distance_man(Point p){
		double dis = Math.abs(p.x - x) + Math.abs(p.y - y);
		return dis;
	}
}
