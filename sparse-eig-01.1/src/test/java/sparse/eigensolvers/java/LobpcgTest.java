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

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import org.junit.Test;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.QR;
import no.uib.cipr.matrix.SVD;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.io.MatrixVectorReader;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import sparse.eigenvolvers.java.*;

public class LobpcgTest extends TestCase  {
	// seed>=0 to generate identical pseudo-random matrices for each run
	// Use seed<0 to generate different pseudo-random matrices for each run
	long seed=10000; 

	// Test dense matrix
	@Test
	public void testDenseMatrix() {
		// Setup problem
		int blockSize=20;
		int n=500;
		int maxIterations=200;
		int verbosityLevel=0;
		DenseMatrix blockVectorX;
		
		// Set A to reproducible matrix for testing purposes
		DenseMatrix A=Utilities.denseMatrixRandomSeed(n,n,seed);
		DenseMatrix Temp=A.copy();
		A.add(Temp.transpose()); //Make it symmetric
		Lobpcg lobpcg=new Lobpcg();
		
		// Set X to reproducible matrix for testing purposes
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		// Run lobpcg to compute smallest eigenvalues and make sure that we converge
		assertTrue(lobpcg.runLobpcg(blockVectorX,A,maxIterations,verbosityLevel)>0);

		// Set X to reproducible matrix for testing purposes
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		// Now run lobpcg again to compute largest eigenvalues and make sure that we converge
		A.scale(-1);
		assertTrue(lobpcg.runLobpcg(blockVectorX,A,maxIterations,verbosityLevel)>0);
	}
	
	// Test Laplacian against theoretical values
	@Test
	public void testLaplacian(){
		int n;
		int verbosityLevel=0;
		int blockSize;
		int maxIterations;
		double tolerance=1e-8;
		int innerIterations;
		int nx,ny,nz;
		double error1,error2,error3;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		
		// Problem setup no preconditioner, Laplacian 20 x 20
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=200;
		nx=20;ny=20;
		CompRowMatrix ACompRow=Utilities.getLaplacian(nx,ny);
		Operator operA=new Operator(ACompRow);
		n=operA.getOperatorSize();
		
		// Set X to reproducible matrix for testing purposes
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,tolerance,maxIterations,verbosityLevel)>0);
		error1=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny),
			lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		assertTrue(error1<1e-8);
		
		// Problem setup with preconditioner, Laplacian 20 x 20
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=100;
		innerIterations=4;
		nx=20;ny=20;
		ACompRow=Utilities.getLaplacian(nx,ny);
		operA=new Operator(ACompRow);
		n=operA.getOperatorSize();
		OperatorPrecCG operT= new OperatorPrecCG(ACompRow);
		operT.setCGNumberIterations(innerIterations);
		n=operA.getOperatorSize();
		
		// Set X to reproducible matrix for testing purposes
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel)>0);
		error2=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny),
				lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		assertTrue(error2<1e-8);

		// Problem setup with CG preconditioner, Laplacian 20 x 20 x 20
		lobpcg=new Lobpcg();
		blockSize=20; maxIterations=200;
		innerIterations=4;
		nx=20;ny=20;nz=20;
		ACompRow=Utilities.getLaplacian(nx,ny,nz);
		operA=new Operator(ACompRow);
		n=operA.getOperatorSize();
		operT= new OperatorPrecCG(ACompRow);
		operT.setCGNumberIterations(innerIterations);
		n=operA.getOperatorSize();
		
		// Set X to reproducible matrix for testing purposes
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel)>0);
		error3=Utilities.sub(Utilities.getLaplacianEigenvalues(blockSize,nx,ny,nz),
			lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		assertTrue(error3<1e-8);
	}
	
	// Test lobpcg constraints
	@Test
	public void testConstraints(){
		int nn;
		DenseMatrix Y;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		DenseMatrix eigVal1,eigVal2,Temp;
		Operator operA,operB;

		// Problem setup
		int blockSize=20;
		int maxIterations=300;
		double tolerance=1e-8;
		int verbosityLevel=0;
		int sizeLaplacian=20;
		operA=new Operator(Utilities.getLaplacian(sizeLaplacian,sizeLaplacian));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		int n=A.numRows();
		DenseMatrix B=Utilities.mult(A,A); // Make B=A^2
		operA=new Operator(A);
		operB=new Operator(B);
		OperatorPrecDenseCholesky operT = new OperatorPrecDenseCholesky(A);
		lobpcg=new Lobpcg();
	
		// Test eigenvalue problem with random constraints
		int yColumns=100;
		Y=Utilities.denseMatrixRandomSeed(n,yColumns,seed);
		// Orthonormalize Y
		QR qr=QR.factorize(Y);
		Y=qr.getQ();
		// Setup and execute SVD
		DenseMatrix U;
		try {
			SVD svd=SVD.factorize(Utilities.mult(Y,Utilities.trans(Y)));
			U=svd.getU();
			DenseMatrix Ycomp=Utilities.getMatrix(U, yColumns, n-1);
			// Get singular vectors for orthogonal complement of Y (null space of Y*Y^T)
			// Compute RITZ values for this subspace and store first blockSize of these
			nn=Ycomp.numColumns();
			double[] eigValues=new double[nn];
			Temp=Utilities.mult(Utilities.trans(Ycomp),Utilities.mult(A,Ycomp));
			Utilities.eig(Temp,eigValues);
			eigVal1=new DenseMatrix(blockSize,1);
			for (int i=0; i<blockSize; ++i) eigVal1.set(i,0,eigValues[i]);
			// Compute Ritz values using constraints
			// This operator has Ritz values plus yColumns zeros
			Temp=Utilities.mult(Ycomp,Utilities.trans(Ycomp));
			operA=new Operator(Utilities.mult(Temp,Utilities.mult(A,Temp)));
			blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
			
			// The constraints should remove the zeros and only keep the Ritz values
			assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,null,Y,tolerance,maxIterations,verbosityLevel)>0);
			eigVal2=lobpcg.getEigenvaluesMatrix();
			assertTrue((Utilities.sub(eigVal2,eigVal1)).norm(Norm.Frobenius)<1e-6);
			
		} catch (NotConvergedException e){}
		
		// Test eigenvalue problem with a preconditioner and constraints
		operA=new Operator(A); // Restore operator A
		blockVectorX=Utilities.denseMatrixRandomSeed(n,2*blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel)>0);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		// Now solve constrained problem
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,operT,Y,tolerance,maxIterations,verbosityLevel)>0);
		eigVal2=lobpcg.getEigenvaluesMatrix();
		
		// Test generalized eigenvalue problem with constraints 
		blockVectorX=Utilities.denseMatrixRandomSeed(n,2*blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,operB,tolerance,maxIterations,verbosityLevel)>0);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,n-1,0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		Temp=Utilities.mult(Utilities.trans(Y),Utilities.mult(B,Y));
		// Now solve constrained problem
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,operB,null,Y,tolerance,maxIterations,verbosityLevel)>0);
		eigVal2=lobpcg.getEigenvaluesMatrix();

		// Test generalized eigenvalue problem with preconditioner and constraints 
		blockVectorX=Utilities.denseMatrixRandomSeed(n,2*blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,operB,operT,tolerance,maxIterations,verbosityLevel)>0);
		Y=Utilities.getMatrix(lobpcg.getEigenvectors(),0,n-1,0,blockSize-1);
		eigVal1=Utilities.getMatrix(lobpcg.getEigenvaluesMatrix(),blockSize,2*blockSize-1,0,0);
		nn=Y.numColumns();
		Temp=Utilities.mult(Utilities.trans(Y),Utilities.mult(B,Y));
		// Now solve constrained problem
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,operB,operT,Y,tolerance,maxIterations,verbosityLevel)>0);
		eigVal2=lobpcg.getEigenvaluesMatrix();
	}
	
	// Test generalized eigenvalues against Matlab computed eigenvalues
	@Test
	public void testGeneralizeEigenvalues(){
		// Problem setup
		int blockSize=30;
		int maxIterations=100;
		double tolerance=1e-10;
		int verbosityLevel=0;
		int nx=19,ny=19;
		Operator operA=new Operator(Utilities.getLaplacian(nx,ny));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		DenseMatrix B=Utilities.mult(A,A);  // Make B=A^2 for testing purposes
		Lobpcg lobpcg=new Lobpcg();
		int n=A.numColumns();
		DenseMatrix blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		
		// Matlab results, A=getLaplacian(19,19), B=A^2, eigs(A,B,361)
		double[] eigsMatlab = new double[]{
			1.257742448321386e-001,
		    1.269439931163143e-001,
		    1.269439931163145e-001,
		    1.281357038671146e-001,
		    1.289093412683051e-001,
		    1.289093412683053e-001,
		    1.301384164527073e-001,
		    1.301384164527079e-001,
		    1.316931275637011e-001,
		    1.316931275637013e-001,
		    1.322047263201810e-001,
		    1.329761236266076e-001,
		    1.329761236266079e-001,
		    1.351342761718838e-001,
		    1.351342761718840e-001,
		    1.353255007435682e-001,
		    1.353255007435682e-001,
		    1.366806124658793e-001,
		    1.366806124658796e-001,
		    1.381966011250104e-001,
		    1.389617162008466e-001,
		    1.389617162008471e-001,
		    1.398416145471212e-001,
		    1.398416145471213e-001,
		    1.412891654094929e-001,
		    1.412891654094931e-001,
		    1.422020474578299e-001,
		    1.422020474578300e-001,
		    1.437280619611563e-001,
		    1.437280619611568e-001,
		};
		// Store Matlab computed eigenvalues
		DenseMatrix eigsMatlab2=new DenseMatrix(eigsMatlab.length,1);
		for (int i=0;i<eigsMatlab.length;++i) eigsMatlab2.set(i,0,eigsMatlab[i]);
		// Run lobpcg
		assertTrue(lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel)>0);
		double error=Utilities.sub(eigsMatlab2,lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
		assertTrue(error<1e-10);
	}

	// Solve standard eigenvalue problem. Test Laplacian against MATLAB computed 
	// eigenvalues and eigenvectors for 3D 20 x 20 x 20 Laplacian with zero boundary conditions.
	@Test
	public void testStandardEigenvaluesMM(){
		boolean flag=true;
		int n;
		int blockSize;
		int verbosityLevel=0;
		int maxIterations=500;
		double tolerance=1e-8;
		int nx=20,ny=20,nz=20;
		double error1,error2;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg=new Lobpcg();;
		String file;
		
		DenseMatrix Eig=new DenseMatrix(0,0);
		DenseMatrix EigVec=new DenseMatrix(0,0);
		
		// Get path
		String path=null;
		File dir = new File(".");
        try {
        	path=dir.getCanonicalPath();
        } catch (Exception e) {}
	        
		// Get eigenvalues
        file=path+"\\src\\test\\resources\\Laplacian20x20x20Eig.mtx";
        if (!System.getProperty("os.name").startsWith("Windows")) file = file.replace("\\", "/"); 
		try {
			 MatrixVectorReader reader = new MatrixVectorReader(new FileReader(file));  
			 CompRowMatrix Eig1 = new CompRowMatrix(reader);
			 Eig=Utilities.CompRowMatrixtoDenseMatrix(Eig1);
		} catch (IOException e) {
			fail("Can't find matrix market file: "+file);
			flag=false;
		}

		// Get eigenvectors
		file=path+"\\src\\test\\resources\\Laplacian20x20x20Vec.mtx";
		if (!System.getProperty("os.name").startsWith("Windows")) file = file.replace("\\", "/"); 
		try {
			 MatrixVectorReader reader = new MatrixVectorReader(new FileReader(file));  
			 CompRowMatrix EigVec1 = new CompRowMatrix(reader);
			 EigVec=Utilities.CompRowMatrixtoDenseMatrix(EigVec1);
		} catch (IOException e) {
			fail("Can't find matrix market file: "+file);
			flag=false;
		}
		
		if (flag){
			CompRowMatrix ACompRow=Utilities.getLaplacian(nx,ny,nz);
			Operator operA=new Operator(ACompRow);
			n=operA.getOperatorSize();
			blockSize=EigVec.numColumns();
			
			// Problem setup with CG preconditioner
			OperatorPrecCG operT= new OperatorPrecCG(ACompRow);
			int innerIterations=4;
			operT.setCGNumberIterations(innerIterations);
		
			// Set blockVectorX to reproducible matrix for testing purposes
			blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
			
			assertTrue(lobpcg.runLobpcg(blockVectorX,operA,null,operT,tolerance,maxIterations,verbosityLevel)>0);
			error1=Utilities.sub(Eig,lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
			error2=Subspace.getPrincipleAnglesLargest(EigVec,lobpcg.getEigenvectors());
			assertTrue(error1<1e-8);
			assertTrue(error2<1e-6);
		} else return;
	}	
	
	// Solve generalized eigenvalue problem. Test Laplacian against MATLAB computed 
	// eigenvalues and eigenvectors for 3D 20 x 20 x 20 Laplacian with zero boundary conditions
	// and 2D 64 x 125 Laplacian on the right hand side
	@Test
	public void testGeneralizedEigenvaluesMM(){
		boolean flag=true;
		int n;
		int blockSize;
		int verbosityLevel=0;
		int maxIterations=500;
		double tolerance=1e-8;
		int nx=20,ny=20,nz=20; // 20x20x20=8000
		double error1,error2;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg=new Lobpcg();;
		String file;
		DenseMatrix Eig=new DenseMatrix(0,0);
		DenseMatrix EigVec=new DenseMatrix(0,0);
		
		// Get path
		String path=null;
		File dir = new File(".");
        try {
        	path=dir.getCanonicalPath();
        } catch (Exception e) {}
        
		// Get eigenvalues
		file=path+"\\src\\test\\resources\\LaplacianAB20x20x20Eig.mtx";
		if (!System.getProperty("os.name").startsWith("Windows")) file = file.replace("\\", "/");
		try {
			 MatrixVectorReader reader = new MatrixVectorReader(new FileReader(file));  
			 CompRowMatrix Eig1 = new CompRowMatrix(reader);
			 Eig=Utilities.CompRowMatrixtoDenseMatrix(Eig1);
		} catch (IOException e) {
			fail("Can't find matrix market file: "+file);
			flag=false;
		}

		// Get eigenvectors
		file=path+"\\src\\test\\resources\\LaplacianAB20x20x20Vec.mtx";
		if (!System.getProperty("os.name").startsWith("Windows")) file = file.replace("\\", "/");
		try {
			 MatrixVectorReader reader = new MatrixVectorReader(new FileReader(file));  
			 CompRowMatrix EigVec1 = new CompRowMatrix(reader);
			 EigVec=Utilities.CompRowMatrixtoDenseMatrix(EigVec1);
		} catch (IOException e) {
			fail("Can't find matrix market file: "+file);
			flag=false;
		}
		
		if (flag){
			CompRowMatrix A=Utilities.getLaplacian(nx,ny,nz);
			nx=64;ny=125; // 64 x 125 = 8000
			CompRowMatrix B=Utilities.getLaplacian(nx,ny);
			n=A.numRows();
			blockSize=EigVec.numColumns();
		
			// Set blockVectorX to reproducible matrix for testing purposes
			blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
			
			assertTrue(lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel)>0);
			error1=Utilities.sub(Eig,lobpcg.getEigenvaluesMatrix()).norm(Norm.Frobenius);
			error2=Subspace.getPrincipleAnglesLargest(EigVec,lobpcg.getEigenvectors());
			assertTrue(error1<1e-8);
			assertTrue(error2<1e-5);
		} else return;
	}
	
	// Test various combinations of input parameters
	@Test
	public void testInputParameters(){
		int blockSize=20;
		int maxIterations=200;
		double tolerance=1e-8;
		int verbosityLevel=0;
		int nx=20,ny=20;
		DenseMatrix blockVectorX;
		Lobpcg lobpcg;
		Operator operA;
		
		// Problem setup
		operA=new Operator(Utilities.getLaplacian(nx,ny));
		DenseMatrix A=operA.getDenseOperatorMatrix();
		int n=A.numColumns();
		DenseMatrix B=Utilities.mult(A,A);  // Make B=A^2
		lobpcg=new Lobpcg();
		operA=new Operator(A);
		
		lobpcg.setVerbosityLevel(0);
		lobpcg.setMaxIterations(200);
		assertTrue(lobpcg.runLobpcg(A)>0);
		assertTrue(lobpcg.runLobpcg(5,A)>0);
		assertTrue(lobpcg.runLobpcg(25,A,B)>0);
		assertTrue(lobpcg.runLobpcg(A,B)>0);
		assertTrue(lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel)>0);
		assertTrue(lobpcg.runLobpcg(A,B,tolerance,maxIterations,verbosityLevel)>0);
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,A,B,tolerance,maxIterations,verbosityLevel)>0);
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,operA,B,tolerance,maxIterations,verbosityLevel)>0);
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		assertTrue(lobpcg.runLobpcg(blockVectorX,tolerance,maxIterations,A,B,verbosityLevel)>0);
		
		// Run lobpcg() with no in line parameters
		lobpcg=new Lobpcg();
		blockVectorX=Utilities.denseMatrixRandomSeed(n,blockSize,seed);
		lobpcg.setBlockVectorX(blockVectorX);
		lobpcg.setA(A);
		OperatorPrecDenseCholesky operT = new OperatorPrecDenseCholesky(A);
		lobpcg.setT(operT);
		lobpcg.setMaxIterations(50);
		lobpcg.setResidualTolerance(1e-8);
		lobpcg.setVerbosityLevel(0);
		assertTrue(lobpcg.runLobpcg()>0); 
		
		// run with simplest parameters
		lobpcg=new Lobpcg();
		lobpcg.setVerbosityLevel(0);
		lobpcg.setA(A);
		lobpcg.setMaxIterations(100);
		assertTrue(lobpcg.runLobpcg()>0);
		
		// run with simplest parameters and change block size
		lobpcg=new Lobpcg();
		lobpcg.setVerbosityLevel(0);
		lobpcg.setBlockSize(25);
		lobpcg.setMaxIterations(100);
		lobpcg.setA(A);
		assertTrue(lobpcg.runLobpcg()>0); 
	}

}
