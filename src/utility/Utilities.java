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

}
