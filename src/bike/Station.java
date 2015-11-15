package bike;

import Geo.Point;
import Geo.Point_latlon;

public class Station {
	public int id;
	public Point_latlon point_latlon;
	public Point point_xy;
	public String address;
	
	
	public Station(){
		point_latlon = new Point_latlon();
		point_xy = new Point();
	}
	
	public double distance_euc(Station s){
		return point_xy.distance_euc(s.point_xy);
	}
	
	public double distance_man(Station s){
		return point_xy.distance_man(s.point_xy);
	}
}
