/*
 * Copyright (C) 2012 Rico Argentati
 * 
 * This file is part of SPARSE-EIG.
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

package sparse.eigenvolvers.java;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;
import java.io.FileReader;
import java.io.IOException;
import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.io.MatrixVectorReader;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import org.netlib.lapack.LAPACK;

/**
SPARSE-EIG Java Utilities class.
<p>
SPARSE-EIG uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
The Utilities class consists of a set of useful methods for running and testing eigenvalue
problems and implementing numerical linear algebra algorithms.
@author Rico Argentati
*/
public class Utilities {

	/* ------------------------
	   Class variables**
	 * ------------------------ */
	static String printString=SparseEigensolverConstants.PRINT_STRING; // Default for print utilities
	
	/* ------------------------
	   Public Methods
	 * ------------------------ */
	
	/**
     * Set print string for convenient display and debugging
     * 
     * @param s
     *            
     */
	public static void setPrintString(String s){
		printString=s;
	}
	
	/**
     * Print dense matrix
     * 
     * @param A
     *            
     */
	public static void print(DenseMatrix A){
		System.out.printf("\nMatrix size: %d x %d\n",A.numRows(),A.numColumns());
		for(int i=0;i<A.numRows();i++) {
			for(int j=0;j<A.numColumns();j++) {
				System.out.printf(printString+" ",A.get(i,j));
			}
			System.out.printf("\n");
		}
	}
	
	/**
     * Print compressed row matrix
     * 
     * @param A
     *            
     */
	public static void print(CompRowMatrix A){
		System.out.printf("\nMatrix size: %d x %d\n",A.numRows(),A.numColumns());
		for(int i=0;i<A.numRows();i++) {
			for(int j=0;j<A.numColumns();j++) {
				System.out.printf(printString+" ",A.get(i,j));
			}
			System.out.printf("\n");
		}
	}

	/**
     * Print dense vector
     * 
     * @param v
     *            
     */
	
	public static void print(DenseVector v){
		System.out.printf("\nVector size: %d\n",v.size());
		for(int i=0;i<v.size();i++) {
			System.out.printf(printString+"\n",v.get(i));
		}
		System.out.printf("\n");
	}
	
	/**
     * Print 1D double array
     * 
     * @param A
     *            
     */
	public static void print(double[] A){
		System.out.printf("\nArray size: %d x 1\n",A.length);
		for(int i=0;i<A.length;i++) {
			System.out.printf(printString+"\n",A[i]);
		}
		System.out.printf("\n");
	}
	
	/**
     * Print 2D double array
     * 
     * @param A
     *            
     */
	public static void print(double[][] A){
		System.out.printf("\nArray size: %d x %d\n",A.length,A[0].length);
		for(int i=0;i<A.length;i++) {
			for(int j=0;j<A[i].length;j++) {
				System.out.printf(printString+" ",A[i][j]);
			}
			System.out.printf("\n");
		}
		System.out.printf("\n");
	}
	
	/**
	 * Print 2D integer array
	 * 
	 * @param A
	 *            
	 */
	public static void print(int[][] A){
		System.out.printf("\nArray size: %d\n",A.length);
		for(int i=0;i<A.length;i++) {
			for(int j=0;j<A[i].length;j++) System.out.printf("%8d ",A[i][j]);
			System.out.printf("\n");
		}
		System.out.printf("\n");
	}
	
	/**
	 * Print 1D long array
	 * 
	 * @param A
	 *            	
	 */
	public static void print(long[] A){
		System.out.printf("\nArray size: %d x 1\n",A.length);
		for(int i=0;i<A.length;i++) {
			System.out.printf("%d\n",A[i]);
		}
		System.out.printf("\n");
	}
	
	/**
     * Print 1D integer array
     * 
     * @param A
     *            
     */
	public static void print(int[] A){
		System.out.printf("\nArray size: %d x 1\n",A.length);
		for(int i=0;i<A.length;i++) {
			System.out.printf("%d\n",A[i]);
		}
		System.out.printf("\n");
	}
	
	/**
     * Print 1D Integer array
     * 
     * @param A
     *            
     */
	public static void print(Integer[] A){
		System.out.printf("\nMatrix size: %d x 1\n",A.length);
		for(int i=0;i<A.length;i++) {
			System.out.printf("%d\n",A[i]);
		}
		System.out.printf("\n");
	}
	
	/**
     * Print dense matrix size and norm
     * 
     * @param A
     *            
     */
	public static void printMatrixInformation(DenseMatrix A){
		System.out.printf("Matrix size = %d x %d, norm=%e\n",
				A.numRows(),A.numColumns(),A.norm(Norm.Frobenius));
	}
	
	/**
     * Print compressed row matrix size and norm
     * 
     * @param A
     *            
     */
	public static void printMatrixInformation(CompRowMatrix A){
		System.out.printf("Matrix size = %d x %d, norm=%e\n",
				A.numRows(),A.numColumns(),A.norm(Norm.Frobenius));
	}

	/**
     * Returns 1D or 2D or 3D CompRowMatrix Laplacian matrix 
     * with simple Dirichlet boundary conditions. Used for testing purposes.
     * See
     * <a href="http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d">
     * Matlab File - laplacian.m </a>
     * @param nx
     * @param ny
     * @param ny
     * <br>
     * <pre><b>Parameters:</b>  
     *       getLaplacian(nx), OR  getLaplacian(nx,ny),  OR  getLaplacian(nx,ny,nz)</pre>
     * 
     * @return CompRowMatrix
     *            
     */
	 public static CompRowMatrix getLaplacian(int ...args){
		 // Parameters: 1D getLaplacian(nx), 2D getLaplacian(nx,ny), 3D getLaplacian(nx,ny,nz)
		 CompRowMatrix Ix,Iy,Iz,D1x,D1y,D1z;
		 CompRowMatrix A=null;
			 
		 // check parameters
		 if (args.length<1 || args.length>3){
			 System.out.println("getLaplacian: wrong number of parameters "+args.length);
			 System.exit(0);
		 }
		 // handle 1D case
		 int nx=args[0];
		 int[][] nz;
		 if (args.length==1){
			 // handle special case n=1
			 if (nx==1){
				 nz=new int[nx][1]; 
				 A=new CompRowMatrix(nx,nx,nz);
				 A.set(0,0,2);
				 return A;
			 }
			 nz=new int[nx][];
			 for (int i=0;i<nx;i++){
				if (i==0){
					 nz[0]=new int[2];
					 nz[0][0]=i;
					 nz[0][1]=i+1;
				}
				else if (i==(nx-1)){
					nz[i]=new int[2];
					nz[i][0]=i-1;
					nz[i][1]=i;
				}
				else {
					nz[i]=new int[3];
					nz[i][0]=i-1;
					nz[i][1]=i;
					nz[i][2]=i+1;
				}
			 }
			 // build and load A
			 A=new CompRowMatrix(nx,nx,nz);
			 for (int i=0;i<nx;i++){
			 	A.set(i,i,2D);
				if (i==0) A.set(i,i+1,-1D);
				else if (i==(nx-1))A.set(i,i-1,-1D);
				else {
					A.set(i,i-1,-1D);
					A.set(i,i+1,-1D);
				}
			 }
			 return A;
		 }

		 // handle 2D and 3D cases
		 // Build sparse identity matrix Ix
		 nz=new int[nx][];
		 for (int i=0;i<nx;i++){
			 nz[i]=new int[1];
			 nz[i][0]=i;
		 }
		 Ix=new CompRowMatrix(nx,nx,nz);
		 for (int i=0;i<nx;i++) Ix.set(i,i,1D);
			 
		 if (args.length==2){
			 // Build sparse identity matrix Iy
			 int ny=args[1];
			 nz=new int[ny][];
			 for (int i=0;i<ny;i++){
				 nz[i]=new int[1];
				 nz[i][0]=i;
			 }
			 Iy=new CompRowMatrix(ny,ny,nz);
			 for (int i=0;i<ny;i++) Iy.set(i,i,1D);
			 D1x=getLaplacian(nx);
			 D1y=getLaplacian(ny);
			 A = add(kron(Iy,D1x),kron(D1y,Ix));
		 } else if (args.length==3){
			// Build sparse identity matrix Iy and Iz
			 int ny=args[1];
			 nz=new int[ny][];
			 for (int i=0;i<ny;i++){
				 nz[i]=new int[1];
				 nz[i][0]=i;
			 }
			 Iy=new CompRowMatrix(ny,ny,nz);
			 for (int i=0;i<ny;i++) Iy.set(i,i,1D);
			 int nnz=args[2];
			 nz=new int[nnz][];
			 for (int i=0;i<nnz;i++){
				 nz[i]=new int[1];
				 nz[i][0]=i;
		 }
		 Iz=new CompRowMatrix(nnz,nnz,nz);
			 for (int i=0;i<nnz;i++) Iz.set(i,i,1D);
			 D1x=getLaplacian(nx);
			 D1y=getLaplacian(ny);
			 D1z=getLaplacian(nnz);
			 A = add(add(kron(Iz, kron(Iy, D1x)),kron(Iz, kron(D1y, Ix))),
				 kron(kron(D1z,Iy),Ix));
		 }
		 return A;
	 }
	 
	 /**
	  * Returns <code>count</code> exact eigenvalues for  1D or 2D or 3D Laplacian matrix
	  * with simple Dirichlet boundary conditions. Used for testing purposes.
	  * See
	  * <a href="http://www.mathworks.com/matlabcentral/fileexchange/27279-laplacian-in-1d-2d-or-3d">
	  * Matlab File - laplacian.m </a>
	  * <p>
	  * See Iterative Methods for Solving Linear Systems, Anne Greenbaum, SIAM 1997, p. 132
	  * <br>
	  * <pre><b>Parameters:</b>  
	  * getLaplacian(count,nx), OR  getLaplacian(count,nx,ny),  OR  getLaplacian(count,nx,ny,nz)</pre>
	  * 
	  * @return DenseMatrix count x 1
	  *            
	  */
	 public static DenseMatrix getLaplacianEigenvalues(int count,int ...args){
		 // getLaplacianEigenvalues(count,nx,ny,nz)
		 int nx=0,ny=0,nz=0,index=0;
		 double ax=0,ay=0,az=0;
		 double[] eigs=new double[0];
		 // check parameters
		 if (args.length<1 || args.length>3){
			 System.out.println("getLaplacianEigenvalues: wrong number of parameters "+args.length);
			 System.exit(0);
		 }
		 // handle 1D case
		 if (args.length==1) {
			 nx=args[0];
			 if (count>nx) {
				 System.out.println("getLaplacianEigenvalues: too many eigenvalues "+count+">"+nx);
				 System.exit(0);
			 }
			 eigs= new double[nx];
			 for (int i=0;i<nx;i++){
				eigs[i]=4*Math.pow(Math.sin((i+1)*Math.PI/(2*(nx+1))),2); 
			 }
		 }
		// handle 2D case
		 if (args.length==2){
			 nx=args[0];
			 ny=args[1];
			 eigs= new double[nx*ny];
			 if (count>nx*ny) {
				 System.out.println("getLaplacianEigenvalues: too many eigenvalues "+count+">"+nx*ny);
				 System.exit(0);
			 }
			 index=0;
			 for (int i=0;i<nx;i++) {
				 ax=4*Math.pow(Math.sin((i+1)*Math.PI/(2*(nx+1))),2); 
				 for (int j=0;j<ny;j++) {
					 ay=4*Math.pow(Math.sin((j+1)*Math.PI/(2*(ny+1))),2);
					 eigs[index++]=ax+ay;
				 }
			 }
		 }
		 // handle 3D case
		 if (args.length==3){
			 nx=args[0];
			 ny=args[1];
			 nz=args[2];
			 eigs= new double[nx*ny*nz];
			 if (count>nx*ny*nz) {
				 System.out.println("getLaplacianEigenvalues: too many eigenvalues "+count+">"+nx*ny*nz);
				 System.exit(0);
			 }
			 index=0;
			 for (int i=0;i<nx;i++) {
				 ax=4*Math.pow(Math.sin((i+1)*Math.PI/(2*(nx+1))),2); 
				 for (int j=0;j<ny;j++) {
					 ay=4*Math.pow(Math.sin((j+1)*Math.PI/(2*(ny+1))),2);
					 for (int k=0;k<nz;k++) {
						 az=4*Math.pow(Math.sin((k+1)*Math.PI/(2*(nz+1))),2);
						 eigs[index++]=ax+ay+az;
					 }
				 }
			 }
		 }

		 // get the requested number of eigenvalues
		 // and store in count x 1 matrix
		 Arrays.sort(eigs);
		 DenseMatrix Temp=new DenseMatrix(count,1);
		 for (int i=0;i<count;++i) Temp.set(i,0,eigs[i]);
		 return Temp;
	 }
	 
	 /**
	  * Add two sparse compressed row format matrices
	  * @param A CompRowMatrix
	  * @param B CompRowMatrix
	  * 
	  * @return compressed row format matrix C=A+B
	  *            
	  */
	 public static CompRowMatrix add(CompRowMatrix A,CompRowMatrix B){
		int[][] nzA=getNz(A);
		int[][] nzB=getNz(B);
		int[][] nzC=new int[nzA.length][];
		int k;
		TreeSet <Integer> t;
		
		// combine nzA and nzB to get nzC 
		for (int i=0;i<nzA.length;i++){
			t = new TreeSet<Integer>();
			k=0; 
			while (k<nzA[i].length) t.add(nzA[i][k++]);
			k=0; 
			while (k<nzB[i].length)	t.add(nzB[i][k++]);
			nzC[i]=new int[t.size()];
			Iterator it =t.iterator();
			k=0;
			while(it.hasNext()) nzC[i][k++]=((Integer) it.next()).intValue();
		}
		CompRowMatrix C=new CompRowMatrix(A.numRows(),A.numColumns(),nzC);
		for (int i=0;i<A.numRows();i++) {
			if (nzC[i].length>0){
				for (int j=0;j<nzC[i].length;j++)
					C.set(i,nzC[i][j],A.get(i,nzC[i][j])+B.get(i,nzC[i][j]));
			}
		 }
		 return C;
	 }
	 
	 /**
	  * Subtract two sparse compressed row format matrices
	  * @param A CompRowMatrix
	  * @param B CompRowMatrix
	  * 
	  * @return compressed row format matrix C=A-B
	  *            
	  */
	 public static CompRowMatrix sub(CompRowMatrix A,CompRowMatrix B){
		int[][] nzA=getNz(A);
		int[][] nzB=getNz(B);
		int[][] nzC=new int[nzA.length][];
		int k;
		TreeSet <Integer> t;
		
		// combine nzA and nzB to get nzC 
		for (int i=0;i<nzA.length;i++){
			t = new TreeSet<Integer>();
			k=0; 
			while (k<nzA[i].length) t.add(nzA[i][k++]);
			k=0; 
			while (k<nzB[i].length)	t.add(nzB[i][k++]);
			nzC[i]=new int[t.size()];
			Iterator it =t.iterator();
			k=0;
			while(it.hasNext()) nzC[i][k++]=((Integer) it.next()).intValue();
		}
		CompRowMatrix C=new CompRowMatrix(A.numRows(),A.numColumns(),nzC);
		for (int i=0;i<A.numRows();i++) {
			if (nzC[i].length>0){
				for (int j=0;j<nzC[i].length;j++)
					C.set(i,nzC[i][j],A.get(i,nzC[i][j])-B.get(i,nzC[i][j]));
			}
		 }
		 return C;
	 }

	 /**
	  * Multiply compressed row matrix by a constant
	  * @param A CompRowMatrix
	  * @param a constant to multiply A by
	  * 
	  * @return compressed row format matrix A=a*A
	  *            
	  */
	 public static CompRowMatrix mult(CompRowMatrix A,double a){
		 for (MatrixEntry e : A)  e.set(e.get()*a); 
		 return A;
	 }
 
	 
	 /**
	  * Kronecker product of sparse compressed row format matrices
	  * @param A CompRowMatrix
	  * @param B CompRowMatrix
	  * 
	  * @return compressed row format matrix
	  *            
	  */	 
	 // Kronecker product of A times B, C = A kron B (all sparse matrices in compressed row format)
	 public static CompRowMatrix kron(CompRowMatrix A,CompRowMatrix B){
		 int zrow,zcol,tempIdx;
		 double value;
		 // load nz pattern for A and B 
		 int[][] nzA;
		 int[][] nzB;
		 nzA=getNz(A);
		 nzB=getNz(B);
		 
		 // create nz patten for kronecker product C=A kron B
		 int[][] nzC=new int[A.numRows()*B.numRows()][];
		 for (int i=0;i<nzC.length;i++)  nzC[i]=new int[0];
		 zrow=0;
		 for (int i=0;i<A.numRows();i++) {
			 for (int j=0;j<B.numRows();j++) {
				nzC[zrow++]=new int[nzA[i].length*nzB[j].length]; 
			 }
		 }
		 zrow=0;
		 for (int i=0;i<A.numRows();i++) {
			 for (int j=0;j<B.numRows();j++) {
				 if (nzC[zrow].length>0){
					 zcol=0;
					 for (int k=0;k<nzA[i].length;k++){
						 tempIdx=nzA[i][k]*B.numColumns();
						 for (int p=0;p<nzB[j].length;p++){
							 nzC[zrow][zcol++]=tempIdx+nzB[j][p]; 
						 }
					 }
				 }
				 zrow++;
			 }
		 }
		 
		 // create compressed row matrix C and load data
		 CompRowMatrix C=new CompRowMatrix(A.numRows()*B.numRows(),A.numColumns()*B.numColumns(),nzC);
		 zrow=0;
		 for (int i=0;i<A.numRows();i++) {
			 for (int j=0;j<B.numRows();j++) {
				 if (nzC[zrow].length>0){
					 zcol=0;
					 for (int k=0;k<nzA[i].length;k++){
						 value=A.get(i,nzA[i][k]);
						 for (int p=0;p<nzB[j].length;p++){
							 C.set(zrow, nzC[zrow][zcol++], value*B.get(j, nzB[j][p]));
						 }
					 }
				 }
				 zrow++;
			 }
		 }

		 return C;
	 }

	 /**
	  * Read file in matrix market format and load data into sparse
	  * matrix in compressed row format
	  * 
	  * @param file String
	  * @return compressed row format matrix 
	  *            
	  */
	 public static CompRowMatrix readMatrixMarketFile(String file){
		 CompRowMatrix A=new CompRowMatrix(0,0,new int[0][0]);
		 try {
			 MatrixVectorReader reader = new MatrixVectorReader(new FileReader(file));  
			 A = new CompRowMatrix(reader);
		 } catch (IOException e) {
			 // Print out the exception that occurred
			 System.out.println("Unable to read Matrix Market file "+e.getMessage());
		}
		return A;
	}
			 
	/**
	 * Compute machine EPS for float
	 * <br>
	 * Machine epsilon gives an upper bound on the relative error due
	 * to rounding in floating point arithmetic.
	 * See http://en.wikipedia.org/wiki/Machine_epsilon
	 * 
	 * @return machFloatEps 
	 *            
	 */
	public static float calculateMachineEpsilonFloat() {
        float machFloatEps = 1.0f;
        do {
           machFloatEps /= 2.0f;
        }
        while ((float)(1.0 + (machFloatEps/2.0f)) != 1.0);
        // System.out.println( "Calculated Machine epsilon: " + machFloatEps );
        return machFloatEps;
    }
	
	/**
	 * Compute machine EPS for double
	 * <br>
	 * Machine epsilon gives an upper bound on the relative error due
	 * to rounding in floating point arithmetic.
	 * See http://en.wikipedia.org/wiki/Machine_epsilon
	 * 
	 * @return machDoubleEps
	 *            
	 */
	public static double calculateMachineEpsilonDouble() {
        double machDoubleEps = 1.0d;
        do {
           machDoubleEps /= 2.0d;
        }
        while ((double)(1.0 + (machDoubleEps/2.0d)) != 1.0);
        // System.out.println( "Calculated Machine epsilon: " + machDoubleEps );
        return machDoubleEps;
    }	 
	 
	
	/**
	 * Compute eigenvalues for symmetric dense matrix
	 * @param A DenseMatrix
	 * 
	 * @return double array of eigenvalues
	 *            
	 */
	public static double[] getEigenvalues(DenseMatrix A){
		double[] eigValues;
		eigValues=new double[A.numColumns()];
		eig(A,eigValues);
		return(eigValues);
	}
	
	/**
	 * Method to solve a generalized eigenvalue problem A*x=lambda*B*x using LAPACK - DSYGV.
	 * A must be symmetric and B must be symmetric positive definite
	 * 
	 * @param A DenseMatrix
	 * @param B DenseMatrix
	 * @param eigValues double array
	 * 
	 */
	public static int eig(DenseMatrix A, DenseMatrix B, double[] eigValues){
		// Get sizes of input matrices
		int n = A.numColumns();
		// Check sizes of input matrices
		if ( A.numRows()!= n | B.numColumns()!=n | B.numColumns()!=n){
			System.out.println("In computeEigDSYGV A and B are not square or not the same size");
			System.exit(0);	
		}
		// Setup all variables
		int itype = 1;
	    String jobz = new String("V"); //Compute eigenvalues and eigenvectors 
	    String uplo = new String("U");
	    int lwork = n*n+1; // Need +1 to handles cases n=1 or n=2
	    double []work = new double[lwork];
	    org.netlib.util.intW info = new org.netlib.util.intW(0);
		// Compute eigenvalues using dsygv
	    LAPACK.getInstance().dsygv(itype,jobz,uplo,n,A.getData(),n,B.getData(),n,eigValues,work,lwork,info);
	    return (info.val);
	}
	
	/**
	 * Method to solve eigenvalue problem A*x=lambda*x using LAPACK - DSYGV.
	 * A must be symmetric
	 * 
	 * @param A DenseMatrix
	 * @param eigValues double array
	 * 
	 */
	public static int eig(DenseMatrix A, double[] eigValues){
		// Get sizes of input matrices
		int n = A.numColumns();
		DenseMatrix B=(DenseMatrix) Matrices.identity(n);
		return eig(A,B,eigValues);
	}

	/**
	* Transpose dense matrix: returns transpose of A (A is unchanged)
	* 
	* @param A DenseMatrix
	* @return A^T DenseMatrix
	*            
	*/
	public static DenseMatrix trans(DenseMatrix A){
		DenseMatrix Temp=new DenseMatrix(A.numColumns(),A.numRows());
		A.transpose(Temp);
		return Temp;
	}

	/**
	 * Multiply dense matrices: return A*B (A and B unchanged)
	 * 
	 * @param A DenseMatrix
	 * @param B DenseMatrix
	 * @return A*B DenseMatrix
	 *            
	 */
	public static DenseMatrix mult(DenseMatrix A,DenseMatrix B){
		DenseMatrix Temp=new DenseMatrix(A.numRows(),B.numColumns());
		A.mult(B, Temp);
		return Temp;
	}
	
	/**
	 * Generate reproducible random matrix for testing purposes
	 * Always generate same matrix for given seed
	 * 
	 * @param n integer
	 * @param m integer
	 * @param seed long
	 * @return Temp DenseMatrix
	 *            
	 */
	public static DenseMatrix denseMatrixRandomSeed(int m, int n, long seed){
		DenseMatrix Temp=new DenseMatrix(m,n);
		if (seed<0) Temp=(DenseMatrix) Matrices.random(m,n); //Used for testing
		else {
			Random random = new Random(seed);
			for (int i=0; i<m; i++){
				for (int j=0; j<n; j++) Temp.set(i, j, random.nextDouble());
			}
		}
		return Temp;
	}
	
	/**
	 * Subtract dense matrices: return A-B (A and B unchanged)
	 * 
	 * @param A DenseMatrix
	 * @param B DenseMatrix
	 * @return C=A-B DenseMatrix
	 *            
	 */
	public static DenseMatrix sub(DenseMatrix A,DenseMatrix B){
		DenseMatrix Temp=new DenseMatrix(A.numRows(),A.numColumns());
		Temp=A.copy();
		Temp.add(-1D,B);
		return Temp;
	}
	
	/**
	 * Convert comress row matrix to dense matrix
	 * 
	 * @param A DenseMatrix
	 * @return Temp DenseMatrix
	 *            
	 */
	public static DenseMatrix CompRowMatrixtoDenseMatrix(CompRowMatrix A){
		DenseMatrix Temp=new DenseMatrix(A.numRows(),A.numColumns());
		for (int i=0;i<A.numRows();i++) {
			 for (int j=0;j<A.numColumns();j++) {
				 Temp.set(i, j, A.get(i, j));
			 }
		 }
		return Temp;
	}
	
	/**
	 * Store specified columns of A in B according to index idx.
	 * Only those columns are returned in B
	 * 
	 * @param A DenseMatrix
	 * @param idx integer array
	 * @return B DenseMatrix
	 *            
	 */
	public static DenseMatrix getMatrix(DenseMatrix A, int[] idx){
		DenseMatrix B=new DenseMatrix(A.numRows(),idx.length);
		// Store specified columns of A in B
		for (int j=0;j<idx.length;++j){
			for (int i=0;i<A.numRows();++i){
				B.set(i,j,A.get(i,idx[j]));
			}
		}
		return B;
	}
	
	/**
	 * Store specified columns of A in B according to interval  j1 to j2 inclusive
	 * Only those columns are returned in B
	 * 
	 * @param A DenseMatrix
	 * @param j1 integer
	 * @param j2 integer
	 * @return B DenseMatrix
	 *            
	 */
	public static DenseMatrix getMatrix(DenseMatrix A, int j1,int j2){
		DenseMatrix B=new DenseMatrix(A.numRows(),j2-j1+1);
		// Store specified columns of A in B
		for (int i=0;i<A.numRows();++i){
			for (int j=0;j<j2-j1+1;++j){
				B.set(i,j,A.get(i,j1+j));
			}
		}
		return B;
	}
	
	/**
	 * Store specified columns of A in B according to rows i1 to i2 and columns j1 to j2 inclusive
	 * Only those columns are returned in B
	 * 
	 * @param A DenseMatrix
	 * @param i1 integer
	 * @param i2 integer
	 * @param j1 integer
	 * @param j2 integer
	 * @return B DenseMatrix
	 *            
	 */
	public static DenseMatrix getMatrix(DenseMatrix A, int i1,int i2,int j1,int j2){
		DenseMatrix B=new DenseMatrix(i2-i1+1,j2-j1+1);
		// Store specified columns of A in B
		for (int i=0;i<i2-i1+1;++i){
			for (int j=0;j<j2-j1+1;++j){
				B.set(i,j,A.get(i1+i,j1+j));
			}
		}
		return B;
	}
	
	/**
	 * Orthonormalize X
	 * Works similar to MATLAB orth.m
	 * Returns matrix who's orthonormal columns span the column space of X
	 * The returned matrix may have fewer columns than X, if X is not full rank
	 * 
	 * @param X DenseMatrix
	 * @param tolUser double (optional argument)
	 * @return DenseMatrix
	 *            
	 */
	public static DenseMatrix orth(DenseMatrix X, double ... tolUser){
		DenseMatrix Temp=new DenseMatrix(0,0);
		double tol;
		if (X.numColumns()==0) return Temp;
		try {
			SVD svd=SVD.factorize(X);
			DenseMatrix U=svd.getU();
			double[] s=svd.getS();
			if (tolUser.length==0) tol=Math.max(X.numRows(), X.numColumns())*s[0]*calculateMachineEpsilonDouble();
			else tol=tolUser[0];
			int ncol=0;
			for (double d : s) if (d>tol) ++ncol;
			Temp=new DenseMatrix(X.numRows(),ncol);
			for (int i=0;i<X.numRows();++i){
				for (int j=0;j<ncol;++j) Temp.set(i,j,U.get(i,j));
			}
		} catch (NotConvergedException e){}
		return Temp;
	}
	
		
	/* ------------------------
	  Package Private Methods
	* ------------------------ */
	 
	// compute non-zero pattern for sparse compressed row matrix 
	static int[][] getNz(CompRowMatrix A){
		int[] row,col;
		int zrow;
			 
		// load nz pattern for A
		row=A.getRowPointers();
		col=A.getColumnIndices();
		int[][] nz=new int[A.numRows()][];
		zrow=0;
		for (int i=0;i<A.numRows();i++)  nz[i]=new int[0];
		for (int i=0;i<row.length-1;i++){
			nz[i]=new int[row[i+1]-row[i]];
			for (int j=0;j<nz[i].length;j++) nz[i][j]=col[zrow++];
		}
		return nz;
	}
	 
	// Compute gram matrix: X^T*X
	static DenseMatrix gramMatrix(DenseMatrix X){
		return mult(trans(X),X);
	}

	// Compute gram matrix using matrix A:  X^T*A*X
	static DenseMatrix gramMatrix(DenseMatrix A,DenseMatrix X){
		return mult(trans(X),mult(A,X));
	}
	
	// Compute gram matrix using operator A:  X^T*A*X
	static DenseMatrix gramMatrix(Operator operA,DenseMatrix X){
		DenseMatrix AX=operA.operatorAction(X); // AX= A*X
		return mult(trans(X),AX);
	}
	
	// Multiply dense arrays: return A*B (A and B unchanged) using basic algorithm
	// This code has been adopted from JAMA
	// Experimental for testing performance
	static double[][] mult(double[][] A,double[][] B){
		int n=A.length;
		int m=A[0].length;
		int p=B[0].length;
		double[][] C=new double[n][p];
		double[] Bcolj = new double[m];
		double[] Arowi;
		double s;
		
		 for (int j = 0; j < p; j++) {
	         for (int k = 0; k < m; k++) {
	            Bcolj[k] = B[k][j];
	         }
	         for (int i = 0; i < n; i++) {
	            Arowi = A[i];
	            s = 0;
	            for (int k = 0; k < m; k++) {
	               s += Arowi[k]*Bcolj[k];
	            }
	            C[i][j] = s;
	         }
	      }
		return C;
	}
	
	// Multiply dense matrices: return A*B (A and B unchanged) using basic algorithm
	// Experimental for testing performance
	static DenseMatrix mult2(DenseMatrix A,DenseMatrix B){
		int n=A.numRows();
		int m=A.numColumns();
		int p=B.numColumns();
		double sum;
		DenseMatrix C=new DenseMatrix(n,p);
		for (int i=0;i<n;i++) {
			 for (int j=0;j<p;j++) {
				 sum=0;
				 for (int k=0;k<m;k++) {
					 sum=sum+A.get(i,k)*B.get(k,j);
				 }
				 C.set(i,j,sum);
			 }
		}
		return C;
	}
	
	// Multiply dense matrices: return A*B (A and B unchanged) using efficient algorithm
	// Experimental for testing performance
	static DenseMatrix mult3(DenseMatrix A,DenseMatrix B){
		double s;
		int n=A.numRows();
		int m=A.numColumns();
		int p=B.numColumns();
		double[] AA=A.getData();
		double[] BB=B.getData();;
		double[] CC=new double[n*p];
		double[] Bcolj = new double[m];
		double[][] Ctemp=new double[n][p];
		for (int j = 0; j < p; j++) {
			int baseBcolj=j*m;
	        for (int k = 0; k < m; k++) {
	        	Bcolj[k] = BB[baseBcolj+k];
	        }
	        for (int i = 0; i < n; i++) {
	        	s = 0;
	            for (int k = 0; k < m; k++) s += AA[i+k*n]*Bcolj[k];
	            CC[i+j*n] = s;
	        }
	    }
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < p; j++) Ctemp[i][j]=CC[i+j*n];
		}
		DenseMatrix C=new DenseMatrix(Ctemp);
		return C;
	}
	
	// Multiply dense matrices: return A*B (A and B unchanged) using efficient algorithm
	// This algorithm will require more memory
	// Experimental for testing performance
	static DenseMatrix mult4(DenseMatrix A,DenseMatrix B){
		double[][] AA=matrixToDoubleArray(A);
		double[][] BB=matrixToDoubleArray(B);
		double[][] CC=mult(AA,BB);
		DenseMatrix C=new DenseMatrix(CC);
		return C;
	}

	// Convert dense matrix to double array
	static double[][] matrixToDoubleArray(DenseMatrix A){
		int n=A.numRows();
		int m=A.numColumns();
		double[] AA=A.getData();
		double[][] TempArray=new double[n][m];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) TempArray[i][j]=AA[i+j*n];
		}
		return TempArray;
	}
	
	// Add dense matrices: return A+B (A and B unchanged)
	static DenseMatrix add(DenseMatrix A,DenseMatrix B){
		DenseMatrix Temp=A.copy();
		return (DenseMatrix) Temp.add(B);
	}
	
	// Add dense matrices: return A=A+B
	static DenseMatrix addEquals(DenseMatrix A,DenseMatrix B){
		A.add(B);
		return A;
	}
	
	// Subtract dense matrices: return A=A-B
	static DenseMatrix subEquals(DenseMatrix A,DenseMatrix B){
		A.add(-1D,B);
		return A;
	}
	
	// Symmetrize square matrix A (A=(A^T+A)/2)
	static DenseMatrix symmetrizeMatrix(DenseMatrix A){
		return (DenseMatrix) add(trans(A),A).scale(.5D);
	}
	
	// Multiply A by scalar  A=a*A
	static DenseMatrix scale(DenseMatrix A,double a){
		return (DenseMatrix) A.scale(a);
	}

	// Method for building an n x m rectangular block MTJ matrix from n*m input matrices
	// Input is row by row. Use null to get zero block.
	// Example: 3 x 3 block matrix: A=blockMatrix(3,3,M1,M2,M3,M4,M5.M6,M7,M8,M9);
	// Example: 3 x 3 block upper triangular matrix:
	// A=blockMatrix(3,3,M1,M2,M3,null,M4.M5,null,null,M6);
	static DenseMatrix blockMatrix(int n, int m, DenseMatrix... M){
		if (n*m != M.length){
			System.out.println("Problem with blockMatrix parameters: n="+n+" m="+m+" matrix count="+M.length);
			return null;
		}
		// Array points to input matrices row by row
		DenseMatrix[][] AA;
		AA=new DenseMatrix[n][m];
		int index=0;
		for (int i=0;i<n;++i){
			for (int j=0;j<m;++j){
				AA[i][j]=M[index++];
			}
		}
		// Create arrays for sizes of matrices and compute dimension of output matrix nn x mm
		int nn=0;
		int mm=0;
		int[] nDimRow=new int[n];
		for (int i=0;i<n;++i){
			nDimRow[i]=0;
			for (int j=0;j<m;++j){
				if (AA[i][j]!=null){
					if(AA[i][j].numRows() > nDimRow[i]) nDimRow[i]=AA[i][j].numRows();
				}
			}
			nn +=nDimRow[i];
		}
		int[] mDimCol=new int[m];
		for (int j=0;j<m;++j){
			mDimCol[j]=0;
			for (int i=0;i<n;++i){
				if (AA[i][j]!=null){
					if(AA[i][j].numColumns() > mDimCol[j]) mDimCol[j]=AA[i][j].numColumns();
				}
			}
			mm +=mDimCol[j];
		}
		
		// Check dimensions of input matrices to see if they have compatible sizes
		for (int i=0;i<n;++i){
			for (int j=0;j<m;++j){
				if ((AA[i][j]!=null) && (AA[i][j].numRows()!=nDimRow[i])){
					System.out.println("Row dimensions are incompatible");
					return null;
				}
			}
		}
		for (int j=0;j<m;++j){
			for (int i=0;i<n;++i){
				if ((AA[i][j]!=null) && (AA[i][j].numColumns()!=mDimCol[j])){
					System.out.println("Column dimensions are incompatible");
					return null;
				}
			}
		}
		// Store matrices in blocks
		DenseMatrix A=new DenseMatrix(nn,mm);
		int rowIndex=0;
		int colIndex=0;
		for (int i=0;i<n;++i){
			colIndex=0;
			for (int j=0;j<m;++j){
				if (AA[i][j]!=null){
					setSubMatrix(rowIndex,colIndex,A,AA[i][j]);
				}
				colIndex+=mDimCol[j];
			}
			rowIndex+=nDimRow[i];
		}
		return A;
	}
	
	
	// Extract submatrix from dense matrix A
	// I am sure that there is a more efficient way to do this
	static DenseMatrix getSubMatrixTMJ(DenseMatrix A,int row1,int row2,int col1,int col2){
		// Validate parameters
		if (row1<0 || row2<0 || col1<0 || col2<0 || row1>row2 || col1>col2 ||
			row2>=A.numRows()|| col2>=A.numColumns()){
			System.out.println("getBlockMatrixTMJ: input parameters invalid");
			System.out.printf("%d %d %d %d\n",row1,row2,col1,col2);
			System.exit(0);
		}
		// Extract data into new submatrix
		DenseMatrix B=new DenseMatrix(row2-row1+1,col2-col1+1);
		for (int i=0;i<=row2-row1;++i){
			for (int j=0;j<=col2-col1;++j){
				B.set(i,j,A.get(row1+i, col1+j));
			}
		}
		return B;
	}
	
	// Store matrix B in A at offset row and col (0 based)
	// I am sure that there is a more efficient way to do this
	static void setSubMatrix(int row, int col, DenseMatrix A, DenseMatrix B){
		// Validate parameters
		if (row<0 || col<0 || row>=A.numRows()|| col>=A.numColumns()){
			System.out.println("setSubMatrix: input parameters invalid");
			System.out.printf("%d %d\n",row,col);
			System.exit(0);
		}
		// Store B in A
		for (int i=0;i<B.numRows();++i){
			for (int j=0;j<B.numColumns();++j){
				A.set(row+i,col+j,B.get(i,j));
			}
		}
	}
	
	// Store double array in jth column of A
	static DenseMatrix setMatrix(DenseMatrix A, double[] x, int j){
		for (int i=0;i<x.length;++i) A.set(i, j, x[i]);
		return A;
	}
	
	// Store columns of B in specified columns of A according to index idx
	static DenseMatrix setMatrix(DenseMatrix A,DenseMatrix B,int[] idx){
		// Store columns of B in specified columns of A
		for (int j=0;j<idx.length;++j){
			for (int i=0;i<A.numRows();++i){
				A.set(i,idx[j],B.get(i,j));
			}
		}
		return A;
	}
	
	// Given gramXBX=X^T*B*X, B-orthonormalize X
	// If B is identity this still works with gramXBX=X^T*X
	static DenseMatrix setOrthCholesky(DenseMatrix gramXBX,DenseMatrix X){
		DenseMatrix Temp=new DenseMatrix(gramXBX.numColumns(),X.numRows());
		DenseCholesky ch = new DenseCholesky(gramXBX.numColumns(),false);
		LowerSPDDenseMatrix LL= new LowerSPDDenseMatrix(gramXBX);
		ch.factor(LL);
		ch.getL().solve(Utilities.trans(X),Temp); // L*Temp=X^T
		return Utilities.trans(Temp);
	}
	
	// Convert LowerSPDDenseMatrix to a DenseMatrix
	static DenseMatrix convertLowerSPDDenseMatrixToDenseMatrix(LowerSPDDenseMatrix A){
		DenseMatrix B=new DenseMatrix(A.numColumns(),A.numRows());
		for (int i=0;i<A.numRows();++i){
			for (int j=0;j<=i;++j){
				B.set(i,j,A.get(i,j));
			}
		}
		return B;
	}
	
	// Compute condition number of symmetric matrix 
	static double getConditionNumber(DenseMatrix A){
		double[] eig=getEigenvalues(A.copy());
		if (eig[0]!=0D) return  Math.abs(eig[eig.length-1]/eig[0]);
		else return Double.NaN;
	}

	// Check orthonormalization
	static double checkOrth(DenseMatrix X){
		return sub(mult(trans(X),X),
			(DenseMatrix) Matrices.identity(X.numColumns())).norm(Norm.Frobenius);
	}
	
	// Check to see if  X and Y are orthogonal
	static double checkOrth(DenseMatrix X,DenseMatrix Y){
		return mult(trans(X),Y).norm(Norm.Frobenius);
	}
	
	// Compute memory usage
	static long getMemoryUsage(){
		Runtime runtime = Runtime.getRuntime(); 
		long memoryInUse=runtime.totalMemory();
		return memoryInUse;
	}
	
	// Compute percentage memory usage
	static double getPercentMemoryUsage(){
		Runtime runtime = Runtime.getRuntime(); 
		long maxMemory = runtime.maxMemory(); 
		long memoryInUse=runtime.totalMemory();
		return memoryInUse*100./maxMemory;
	}
	
	// Get the rank of a dense matrix
	static int getRank(DenseMatrix A){
		double tol = Math.sqrt(Utilities.calculateMachineEpsilonDouble())*A.numRows();
		QR qr=QR.factorize(A.copy());
		UpperTriangDenseMatrix R=qr.getR();
		int n=Math.min(R.numRows(), R.numColumns());
		int rank=0;
		for (int i=0;i<n;i++) if (Math.abs(R.get(i,i))>tol) ++rank;
		return rank;
	}
	
	// Get a dense vector from column of A
	static DenseVector getVectorMatrix(DenseMatrix A,int col){
		DenseVector v=new DenseVector(A.numRows());
		for (int i=0;i<A.numRows();i++) v.set(i,A.get(i, col));
		return v;
	}
	
	// Set column of A to dense vector
	static void setVectorMatrix(DenseMatrix A,DenseVector v,int col){
		for (int i=0;i<A.numRows();i++) A.set(i,col,v.get(i));
	}
	
}

