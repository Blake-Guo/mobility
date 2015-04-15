/*
 * Copyright (C) 2012 Rico Argentati
 * 
 * This file is part of SEJ (Sparse Eigensolvers for Java).
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package sparse.eigensolvers.java;

//Import MTJ classes
//see http://code.google.com/p/matrix-toolkits-java/
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

// This example illustrates how to find eigenvalues for a symmetric dense matrix
// Loppcg is really designed to handle sparse matrices, but this example is included as an illustration.
// There are many java programs in existence to compute eigenvalues for dense matrices.
public class ExampleDenseMatrix {
	public static void main(String[] args) {
		// Setup problem
		int n=500;
		int maxIterations=200;
		DenseMatrix A=(DenseMatrix) Matrices.random(n,n);
		//Make it symmetric
		DenseMatrix Temp=A.copy();
		A.add(Temp.transpose());
		Lobpcg lobpcg=new Lobpcg();
//		lobpcg.setDebugOn();
		
		// Run lobpcg with mostly default parameters
		lobpcg.runLobpcg(A,maxIterations);
		
		// Now eigenvectors are transferred to blockVectorX
		DenseMatrix blockVectorX=lobpcg.getEigenvectors();
		// Can obtain eigenvalues for further processing in dense matrix
		Utilities.print(lobpcg.getEigenvaluesMatrix());
		// Or get eigenvalues in double array
		Utilities.print(lobpcg.getEigenvalues());
		
		// Now run lobpcg again to compute largest eigenvalues
		int verbosityLevel=0;
		A.scale(-1);
		lobpcg.runLobpcg(A,maxIterations,verbosityLevel); // no output
		// Print negative of largest eigenvalues
		Utilities.print(lobpcg.getEigenvaluesMatrix());
	}
}
