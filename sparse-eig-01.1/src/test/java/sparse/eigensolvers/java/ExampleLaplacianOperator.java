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
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

//Import sparse eigensolver package
//see http://code.google.com/p/sparse-eigensolvers-java/
import sparse.eigenvolvers.java.*;

public class ExampleLaplacianOperator {
	
	// This example illustrates how to find eigenvalues for a sparse operator stored
	// in compressed row format.
	// We generate a 3D Laplacian for testing purposes and check the computed eigenvalues
	// against those computed theoretically.
	public static void main(String[] args) {
		// Setup problem
		int blockSize=20;
		int maxIterations=100;
		double tolerance=1e-8;
		int verbosityLevel=1;
		int innerIterations=5;
		int nx=20,ny=20,nz=20;
		double error;
		Lobpcg lobpcg=new Lobpcg();
		CompRowMatrix ACompRow=Utilities.getLaplacian(nx,ny,nz); //Laplacian 20 x 20 x 20
		int n=ACompRow.numColumns();
		DenseMatrix blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
				
		// Run without preconditioner which takes many more iterations
		maxIterations=300;
		lobpcg.runLobpcg(blockVectorX,ACompRow,tolerance,maxIterations,verbosityLevel);
		
		// Run lobpcg with preconditioner which takes many fewer iterations
		// We need place holder "null" since we are not solving
		// a generalized eigenvalue problem and there is no operB
		// Use conjugate gradient preconditioner (operT approximates inverse of ACompRow)
		OperatorPrecCG operT= new OperatorPrecCG(ACompRow); 
		operT.setCGNumberIterations(innerIterations); // Need to try different values to see best convergence
		maxIterations=80;
		blockVectorX=(DenseMatrix) Matrices.random(n,blockSize);
		lobpcg.runLobpcg(blockVectorX,ACompRow,null,operT,tolerance,maxIterations,verbosityLevel);
		error=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny,nz),
			lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		
		System.out.println("Test results against theoretically computed eigenvalues:");
		System.out.printf("Error= %e (Laplacian 20 x 20 x 20 with preconditioner)\n",error);
	}
}
