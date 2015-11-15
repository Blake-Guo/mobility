package Geo;

public class Point_latlon {
	public double lat;
	public double lon;
	
	public Point_latlon(double latitude, double longitude ){
		this.lat = latitude;
		this.lon = longitude;
	}
	
	
	public Point_latlon(){
		this.lat = 0;
		this.lon = 0;
	}
}
