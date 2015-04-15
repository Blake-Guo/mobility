package sparse.eigensolvers.java;

import java.net.URL;
import java.net.URLEncoder;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import org.junit.Test;


import no.uib.cipr.matrix.io.MatrixVectorReader;
import no.uib.cipr.matrix.sparse.CompRowMatrix;


public class TestMisc {
	public static void main(String[] args) {
		// Get path
		String path=null;
		File dir = new File(".");
        try {
        	path=dir.getCanonicalPath();
        } catch (Exception e) {}
	        
		// Get eigenvalues
        String file1=path+"\\src\\test\\resources\\Laplacian20x20x20Eig.mtx";
        System.out.println(file1);
        if (!System.getProperty("os.name").startsWith("Windows")) file1 = file1.replace("\\", "/"); 
        System.out.println(file1);

	}

	
}
