package spectral;

import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.UpperSymmDenseMatrix;
import no.uib.cipr.matrix.io.MatrixInfo;
import no.uib.cipr.matrix.io.MatrixSize;
import no.uib.cipr.matrix.io.MatrixVectorReader;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.CompRowMatrix;


public class Test  {
	public static void main(String[] args)throws IOException{

//		/*
//		 * 5,0,2
//		 * 0,3,0
//		 * 2,0,1
//		 * */
		int[][] ci = new int[3][];
		ci[0] = new int[]{0,2};
		ci[1] = new int[]{1};
		ci[2] = new int[]{0,2};
		
		CompRowMatrix crsMat = new CompRowMatrix(3,3,ci);
		crsMat.set(0, 0, 5);
		crsMat.set(0, 2, 2);
		crsMat.set(1, 1, 3);
		crsMat.set(2, 0, 2);
		crsMat.set(2, 2, 1);
		
        ArpackSym solver = new ArpackSym(crsMat);
        
        Map<Double, DenseVectorSub> results = solver.solve(2, ArpackSym.Ritz.SA);
        
        for (Map.Entry<Double, DenseVectorSub> e : results.entrySet()) {
        	
        	System.out.println(e.getKey());
        	
        	DenseVectorSub evec = e.getValue();
        	for(int i=0;i<evec.size();i++)
        	{
        		System.out.print(evec.get(i) + " ");
        	}
        	System.out.println();

        }

	}

}
