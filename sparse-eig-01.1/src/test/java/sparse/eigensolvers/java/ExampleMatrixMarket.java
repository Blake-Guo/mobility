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
import no.uib.cipr.matrix.sparse.CompRowMatrix;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

// This example illustrates how user can find several of the largest eigenvalues 
// after loading a matrix market file
// Note that string file needs to be defined based on where user stores matrix market files
// See http://math.nist.gov/MatrixMarket/
public class ExampleMatrixMarket {
	public static void main(String[] args) {
		// Problem setup
		int blockSize=10;
		int maxIterations=100;
		int verbosityLevel=1;
		double tolerance=1e-6;
//		String file="C:\\MatrixMarket\\plat1919.mtx";
		String file="C:\\MatrixMarket\\test1A.mtx";
		Lobpcg lobpcg=new Lobpcg();
		CompRowMatrix A=Utilities.readMatrixMarketFile(file);
//		A.scale(-1); // Multiply operator times -1 to compute the largest eigenvalues
		
		// Run lobpcg
		// Note that if the first parameter is an integer, a random set of starting vectors
		// defining the initial subspace is generated inside lobpcg
		lobpcg.runLobpcg(blockSize,A,tolerance,maxIterations,verbosityLevel);
		
		// Now eigenvectors are transferred to blockVectorX
		DenseMatrix blockVectorX=lobpcg.getEigenvectors();
		// Can obtain eigenvalues for further processing in dense matrix
		Utilities.print(lobpcg.getEigenvaluesMatrix());
		// Or get eigenvalues in double array
		Utilities.print(lobpcg.getEigenvalues());
	}
}
