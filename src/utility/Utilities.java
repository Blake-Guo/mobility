package utility;

import java.io.*;

public class Utilities {
	
	public static int countNumberOfLine(String filePath) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(filePath));
		
		int c = 0;
		String line = "";
		
		while((line = br.readLine()) != null)
		{
			c++;
		}
		
		return c;
	}
	
	public static void norm(double[] a){
		double sum = 0;
		for(int i=0;i<a.length;i++)
			sum += a[i] * a[i];
		
		sum = Math.sqrt(sum);
		
		for(int i=0;i<a.length;i++)
			a[i] /= sum;
	}

}
