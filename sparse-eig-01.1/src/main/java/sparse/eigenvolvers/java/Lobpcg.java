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

package sparse.eigenvolvers.java;

import java.util.ArrayList;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.Matrix.Norm;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
SEJ Java Lobpcg class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
This algorithm is based on BLOPEX (Block Locally Optimal Preconditioned
Eigenvalue Xolvers, see [1]).
This Java implementation is called Lobpcg (Locally Optimal Block Preconditioned
Conjugate Gradient Method) and
computes a set of the smallest eigenvalues lambda and corresponding eigenvectors
X of the generalized eigenproblem A*x=lambda B*x, where 
Hermitian operators A and B are given, as well as an optional preconditioner, T.
The operators B and T must, in addition, be positive definite.
B is optional and if it is not provided the standard eigenvalue problem A*x=lambda*x
is solved.  To compute the largest
eigenpairs of A, simply apply the code to A multiplied by
-1. The code does not involve any matrix factorizations of A and B.
The matrices A, B and T may be provided as dense or compressed row matrices or as
an operator function provided in a completely flexible way by the user.
<p>
If a dense matrix Y is provided then the search is done
in the orthogonal (B-orthogonal if B is given) complement of the column-space of Y.
<p>
The arguments residualTolerance and maxIterations control the residual 
tolerance and maximum number of steps, 
and verbosityLevel = 0, 1, or >1 controls the amount of printed
information. The eigenvalue history of all iterative lambdas, and
the residual norms history  of residuals can be obtained using the 
{@link #getEigenvalueHistory()} and {@link #getResidualNormsHistory()} methods
after {@link #runLobpcg(Object...)} is run.
<p>
This code does not compute all eigenvalues. It is designed to compute no more than ~20%,
and attempting to compute more that 20% of the eigenpairs may cause the code to crash.
<p>
For background see [1], [2] and [3] below. The interface and structure of the code functionality is very
similar to lobpcg.m, the Matlab version [2].
<p>
<DT><B>Example of use:</B></DT>
See {@link #runLobpcg(Object...)} for details concerning parameters. 
<P><PRE>
Example with all possible arguments:
runLobpcg(X,A,B,T,Y,residualTolerance,maxIterations,verbosityLevel)

Example of just one argument (solves A*x=lambda*x):
In this case starting vectors for the initial subspace are generated randomly and the 
default block size is used.
lobpcg.runLobpcg(A);

Example of just two arguments (solves A*x=lambda*x):
lobpcg.runLobpcg(X,A);

Example of X, A, tol, maxit, verb (solves A*x=lambda*x):
lobpcg.runLobpcg(X,A,1e-5,50,1);

Example of X, A, B, tol, maxit, verb (solves A*x=lambda*B*x):
lobpcg.runLobpcg(X,A,B,1e-5,50,1);

Example of no B or T, but there is constraint matrix Y (solves A*x=lambda*x):
lobpcg.runLobpcg(X,A,null,null,Y,1e-5,50,1);

Example of no B but there is preconditioner T (solves A*x=lambda*x):
lobpcg.runLobpcg(X,A,null,T);
</PRE></DD>
</DL>
<DT><B>References:</B></DT>
[1] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov, <a
href="http://www-math.ucdenver.edu/~rargenta/index_files/BlopexSiamSISC29_5.pdf">Block
Locally Optimal Preconditioned Eigenvalue Xolvers (BLOPEX) in hypre and PETSc</a>,
SIAM Journal on Scientific Computing (SISC), Vol. 29 (2007),
No. 5, pp. 2224-2239. 
<br>
[2] Lobpcg solves Hermitian partial generalized eigenproblems using preconditioning - 
<a href="http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m"> Lobpcg (MATLAB)<a>
<br>
[3] C-version of Lobpcg -<a href="http://code.google.com/p/blopex/"> Blopex<a>
<br>
@author Rico Argentati
*/

public class Lobpcg extends SparseEigensolver {
    	
	/* ------------------------
	   Constructors
	 * ------------------------ */
	public Lobpcg(int ... args){
	}
	
	/* ------------------------
	   Public Methods
	 * ------------------------ */
	/**
	 * Method for running Lobpcg with input arguments.
	 * <br>
	 * Example with all possible arguments:<br>
	 * runLobpcg(X,A,B,T,Y,residualTolerance,maxIterations,verbosityLevel)
     *  @param blockVectorX initial approximation to eigenvectors, dense matrix (n-by-blockSize).
     *  X must be full rank. Note that in place of X, the blocksize may be given, then the X is 
     *  generated by lobpcg randomly. 
     *  @param operA the operator of the problem, can be given as a dense matrix, compressed
     *  row matrix or as a user defined operator (function). operA, operB and operT conceptually implement
     *  a matrix vector multiply, but this can be done as a function, where operA, operB 
     *  and operT are n x n operators.
     *  <br>
	 * Optional input:
	 *  @param operB the second operator, can be given as a dense matrix, compressed
     *  row matrix or as an operator
	 *  @param operT the preconditioner, can be given as a dense matrix, compressed
     *  row matrix or as an operator. T should be a approximation to the inverse of A
     *  to solve  A*x=lambda*x
     *  @param blockVectorY  constraints input, is a dense matrix. The iterations will be performed
     *  in the orthogonal (B-orthogonal if B is given) complement of the column-space of Y. 
     *  Y must be full rank and have n rows.
     *  @param residualTolerance default is n*sqrt(eps)
     *  @param maxIterations default is min(n,20)
     *  @param verbosityLevel default is 1 (0 - no output, 1 standard output, >1 verbose output)
     *  
     *  <p> 
     *  Note: if there is a T, but no B use null as place holder for B. If there is a Y, but no
     *  B or T use place holder(s) null.
     *  @return convergence status (1 if converged, -1 if did not converge)
	 *
	 */
	public int runLobpcg(Object ... args){

		// Set defaults
		operT=new Operator();
		operB=new Operator();
		blockVectorY=new DenseMatrix(0,0);

		// Parse command line arguments
		if (args.length<1 || args.length>8){
			System.out.print("Wrong number of arguments: "+ args.length);
			System.exit(0);
		}

		int countMatrixOrOperator=0;
		int countInteger=0;
		for (Object ob:args){
//			System.out.println(args.length+" "+ob.getClass());
			if (ob instanceof Double) setResidualTolerance(((Double) ob).doubleValue());
			else if (countMatrixOrOperator!=0 && ob instanceof Integer) {
				switch (countInteger) {
        			case 0: setMaxIterations(((Integer) ob).intValue());
							++countInteger;
							break;
        			case 1: setVerbosityLevel(((Integer) ob).intValue());
							++countInteger;
							break;
					default: System.out.print("Too many Integer arguments: "+countInteger);
							System.exit(0);
				}
				
			} else {
				// A string or null is used as place holder
				if (ob instanceof String || ob==null){
					++countMatrixOrOperator;
					continue;
				}
				switch (countMatrixOrOperator) {
            		case 0: ++countMatrixOrOperator;
            				saveBlockVectorX=null;
            				setBlockVectorX(new DenseMatrix(0,0));
            				// check for special case where user wants to use default random input matrix
            				if (ob instanceof DenseMatrix) {
            					DenseMatrix Temp=(DenseMatrix) ob;
            					if (Temp.numRows() != Temp.numColumns()){
            						saveBlockVectorX=(DenseMatrix) ob;
            						setBlockVectorX((DenseMatrix) ob);
            						break;
            					}
            				}
            				// if first parameter is integer assume that this is the block size
            				// for random input matrix
            				else if (ob instanceof Integer){
            					setBlockSize((Integer) ob);
            					break;
            				}
            		case 1: if (ob instanceof Operator) setA((Operator) ob);
            				else if (ob instanceof CompRowMatrix) setA((CompRowMatrix) ob);
            				else  setA((DenseMatrix) ob);
            				++countMatrixOrOperator;
            				break;
            		case 2: if (ob instanceof Operator) setB((Operator) ob);
            				else if (ob instanceof CompRowMatrix) setB((CompRowMatrix) ob);
    						else setB((DenseMatrix) ob);
    						++countMatrixOrOperator;
    						break;
            		case 3: if (ob instanceof Operator) setT((Operator) ob);
            				else if (ob instanceof CompRowMatrix) setT((CompRowMatrix) ob);
							else  setT((DenseMatrix) ob);
							++countMatrixOrOperator;
							break;
            		case 4: setBlockVectorY((DenseMatrix) ob);
							break;
					default: System.out.print("Too many Matrix and Operator arguments. Exiting..."
								+countMatrixOrOperator);
							System.exit(0);
				}
			}
		}
		
		// Run lobpcg
		int result=runLobpcg();
		// Update values in input eigenvector matrix
		if (saveBlockVectorX!=null){
			for (int i=0;i<blockVectorX.numRows();i++) {
				for (int j=0;j<blockVectorX.numColumns();j++) {
					saveBlockVectorX.set(i,j,blockVectorX.get(i,j));
				}
			}
		}
		return result;
		
	}

	/**
    * Method for running lobpcg using set methods.
    * Use subset of desired SparseEigensolver set methods to set up parameters.
    * Then run method {@link #runLobpcg()}.
    * <br>
    * {@link #setBlockSize(int blockSize)}<br>
    * {@link #setBlockVectorX(DenseMatrix blockVectorX)}<br>
    * {@link #setA(DenseMatrix A)}<br>
    * {@link #setA(CompRowMatrix A)}<br>
    * {@link #setA(Operator A)}<br>
    * {@link #setB(DenseMatrix B)}<br>
    * {@link #setB(CompRowMatrix B)}<br>
    * {@link #setB(Operator B)}<br>
    * {@link #setT(DenseMatrix T)}<br>
    * {@link #setT(CompRowMatrix T)}<br>
    * {@link #setT(Operator T)}<br>
    * {@link #setBlockVectorY(DenseMatrix blockVectorY)}<br>
    * {@link #setMaxIterations(int maxIterations)}<br>
    * {@link #setResidualTolerance(double residualTolerance)}<br>
    * {@link #setVerbosityLevel(int verbosityLevel)}<br>
	*
	* @return convergence status (1 if converged, -1 if did not converge)
	*/
	public int runLobpcg(){
		// Check for existence of operator A
		// The operator A must exist, we can default on everything else
		if (!operA.getExists()) {
			System.out.print("Operator operA does not exist! Exiting...");
			System.exit(0);
		}
		
		// Check to see if blockVectorX exists, if not create randomly
		if (blockVectorX.numColumns()==0){
			setBlockVectorX((DenseMatrix)Matrices.random(operA.getOperatorSize(),blockSize));
		}
		
		int n=operA.getOperatorSize(); 			// Get problems size
		eigenvalues=new double[blockSize];   		// Vector of eigenvalues
		residualNorms=new double[blockSize];	// Vector of residual norms
		currentBlockSize=blockSize;			    // Block size of non-converged eigenvectors
		int currentSubspaceDimension=0;			// Dimension of current subspace (subspace iteration)
		
		// Check block size
		if (n-blockVectorY.numColumns()<5*blockSize){
			System.out.printf("Block size is too big for problem: %d x %d\n",n,blockSize);
			System.exit(0);
		}

		// Set more defaults
		if (!residualToleranceSet) 
			residualTolerance = Math.sqrt(Utilities.calculateMachineEpsilonDouble())*blockVectorX.numRows();
		if (!maxIterationsSet) maxIterations = Math.min(blockVectorX.numRows(),
			SparseEigensolverConstants.MAX_NUMBER_ITERATIONS);
		
		DenseMatrix blockVectorP=new DenseMatrix(n,blockSize);
		DenseMatrix blockVectorR=new DenseMatrix(n,blockSize);
		DenseMatrix blockVectorAP=new DenseMatrix(n,blockSize);
		DenseMatrix blockVectorAX=new DenseMatrix(n,blockSize);
		DenseMatrix blockVectorBX=new DenseMatrix(0,0);
		DenseMatrix blockVectorBP=new DenseMatrix(0,0);
		DenseMatrix blockVectorAR=new DenseMatrix(n,blockSize);
		DenseMatrix blockVectorBR=new DenseMatrix(0,0);
		DenseMatrix blockVectorBY=new DenseMatrix(0,0);
		DenseMatrix gramA,gramB;
		DenseMatrix gramXAX=null,gramXAR,gramXAP,gramXBX,gramXBR,gramXBP;
		DenseMatrix gramRAR,gramRBP,gramRBR,gramRAP,gramPAP,gramPBP;
		DenseMatrix CP, CX, CR;
		DenseMatrix EIG=new DenseMatrix(blockSize,blockSize); 	// Used to store diagonal matrix of eigenvalues
		int nn; 	// temporary integer
		int iterationCount=0;
		boolean convergenceFlag=false;
		boolean restartFlag=false; // Not used now
		long maxMemoryUsed=0;
		DenseMatrix Temp;
		Integer[] tempInteger;
		int[] mskIdx=null;
		ArrayList<Integer> activeMask=null;
		int rank;
		int restartCount=0;
		double conditionG=0;        // condition number to save
		double conditionGMean=0;    // Log base 10 of condition number plus 1
		double conditionGramA=0;
		double conditionGramB=0;
		int k;

		residualNormsHistory=new DenseMatrix(blockSize,maxIterations);
		eigenvalueHistory=new DenseMatrix(blockSize,maxIterations);
		timeHistory=new long[maxIterations];
		residualMaximumHistory=new double [maxIterations];
		conditionGHistory=new double [maxIterations+1];
		conditionGHistory[0]=-Math.log10(Utilities.calculateMachineEpsilonDouble())/2.0;
		
		if (verbosityLevel>0){
			if (verbosityLevel>1){
				System.out.printf("\nMachine eps: %e \n",Utilities.calculateMachineEpsilonDouble());
				System.out.printf("Number of processors available to the Java Virtual Machine: %d \n",
					Runtime.getRuntime().availableProcessors());
			}
		    System.out.printf("Maximum number of iterations: %d \n",maxIterations);
		    System.out.printf("Tolerance: %e \n",residualTolerance);
		    System.out.printf("Operator size: %d x %d\n",operA.getOperatorSize(),operA.getOperatorSize());
		    System.out.printf("Initial guess has %d colunms and %d rows\n",
		    		blockVectorX.numRows(),blockVectorX.numColumns());
		    if (!operB.getExists()) 
		    	System.out.printf("Solving standard eigenvalue problem, not generalized\n");
		    if (!operT.getExists()) 
		    	System.out.printf("No preconditioner is detected\n");
		    if (blockVectorY.numColumns()==0) 
		    	System.out.printf("No matrix of constraints is detected\n\n");
		}
		
		// Initialize for first iteration or restart
		restartFlag=false;

		// Initialize soft locking mask
		activeMask=new ArrayList<Integer>();
		for (int i=0; i < blockSize; i++) activeMask.add(i);
		tempInteger=new Integer[activeMask.size()];
		tempInteger=activeMask.toArray(tempInteger);
		mskIdx=new int [activeMask.size()];
		for (int i=0; i < tempInteger.length; i++) mskIdx[i]=tempInteger[i].intValue();
		
		// Initialize block matrices for operator B
		if (operB.getExists()) {
			blockVectorBX=new DenseMatrix(n,blockSize);
			blockVectorBP=new DenseMatrix(n,blockSize);
			blockVectorBR=new DenseMatrix(n,blockSize);
		}
		
		// Implement contraints
		// Project X onto subspace B-orthogonal to Y (or orthogonal to Y)
		// P=I-Y*Y^TB is now orthogonal projector onto subspace B-orthogonal to Y
		if (blockVectorY.numColumns() != 0){
			// B orthonormalize Y+iterationNumber
			if (operB.getExists()){
				blockVectorBX=operB.operatorAction(blockVectorX);
				blockVectorBY=operB.operatorAction(blockVectorY);
				Temp=Utilities.mult(Utilities.trans(blockVectorY),blockVectorBY);
				}
			else Temp=Utilities.mult(Utilities.trans(blockVectorY),blockVectorY);
			// Make sure that we are really symmetric 
			Utilities.symmetrizeMatrix(Temp);
			blockVectorY=Utilities.setOrthCholesky(Temp,blockVectorY);
			if (operB.getExists()) blockVectorBY=Utilities.setOrthCholesky(Temp,blockVectorBY);
			
			// project X onto subspace B-orthogonal to Y (or orthogonal to Y)
			// X=X-Y*(Y^T*B*X)
			if (operB.getExists()) blockVectorX=Utilities.sub(blockVectorX,
				Utilities.mult(blockVectorY,Utilities.mult(Utilities.trans(blockVectorBY),blockVectorX)));
			// X=X-Y*(Y^T*X)
			else  blockVectorX=Utilities.sub(blockVectorX,
				Utilities.mult(blockVectorY,Utilities.mult(Utilities.trans(blockVectorY),blockVectorX))); 
		}
			
		// Check to see if full rank
		rank=Utilities.getRank(blockVectorX);
		if (rank<blockVectorX.numColumns()){
			System.out.println("Block vector input not full rank!");
			System.out.printf("Rank = %d, Block size = %d\n",rank,blockVectorX.numColumns());
		}
		
		// B orthonormalize X
		if (operB.getExists()){
			blockVectorBX=operB.operatorAction(blockVectorX);
			gramXBX=Utilities.mult(Utilities.trans(blockVectorX),blockVectorBX);
			}
		else gramXBX=Utilities.mult(Utilities.trans(blockVectorX),blockVectorX);
		Utilities.symmetrizeMatrix(gramXBX); // Make sure that we are really symmetric
		
		blockVectorX=Utilities.setOrthCholesky(gramXBX,blockVectorX);
		if (operB.getExists()) blockVectorBX=Utilities.setOrthCholesky(gramXBX,blockVectorBX);

		// Compute the initial Ritz vectors: solve the eigenproblem
		blockVectorAX=operA.operatorAction(blockVectorX);
		gramXAX=Utilities.mult(Utilities.trans(blockVectorX),blockVectorAX);
		Utilities.symmetrizeMatrix(gramXAX); // Make sure that we are really symmetric 
		
		// Compute eigenvalues and eigenvectors
		eigenvalues=new double[blockSize];
		if (Utilities.eig(gramXAX,eigenvalues)!=0)
			System.out.println("Generalized eigenvalue solution could not be completed (LAPACK - DSYGV)");
			currentSubspaceDimension=gramXAX.numColumns();
		// Assign eigenvalues to diagonal of EIG
		EIG=new DenseMatrix(blockSize,blockSize);
		for (int i=0;i<blockSize;++i){
			EIG.set(i,i,eigenvalues[i]);
		}
		
		// Update X, AX and BX
		blockVectorX=Utilities.mult(blockVectorX,gramXAX);
		blockVectorAX=Utilities.mult(blockVectorAX,gramXAX);
		if (operB.getExists()) blockVectorBX=Utilities.mult(blockVectorBX,gramXAX);
		
		/*****************************************************************
		 *  
		 * Main loop
		 * 
		 ******************************************************************/
		long t0 = System.currentTimeMillis(); // store current time
		long timeIteration;
		for (int iterationNumber=1;iterationNumber<=maxIterations;++iterationNumber){
			iterationCount=iterationNumber;
			// Store current time
			timeIteration=System.currentTimeMillis();
			
			if (iterationNumber==1){
				// Compute residuals
				if (operB.getExists()) blockVectorR=Utilities.sub(blockVectorAX, Utilities.mult(blockVectorBX, EIG));
				else blockVectorR=Utilities.sub(blockVectorAX, Utilities.mult(blockVectorX, EIG));
				// Update active mask and active mask array
				for (int j=0;j<blockSize;++j){
					residualNorms[j]=Utilities.getMatrix(blockVectorR,j,j).norm(Norm.Frobenius);
					if (iterationNumber>1){ 
						if (residualNorms[j]<residualTolerance){
							// Remove index if in list
							if (activeMask.indexOf(j)>-1) activeMask.remove(activeMask.indexOf(j));
						}
					}
				}
				// Update active mask
				tempInteger=new Integer[activeMask.size()];
				tempInteger=activeMask.toArray(tempInteger);
				mskIdx=new int [activeMask.size()];
				for (int i=0; i < tempInteger.length; i++) mskIdx[i]=tempInteger[i].intValue();
			}
			
			currentBlockSize=activeMask.size();
			if (currentBlockSize==0){  // Check for convergence of all eigenvalues
				convergenceFlag=true;
				break; 
			}
			
			// Implement constraints
			// Project R onto subspace B-orthogonal to Y (or orthogonal to Y)
			if (blockVectorY.numColumns() != 0){
				Temp=Utilities.getMatrix(blockVectorR,mskIdx);
				// R=R-Y*(Y^T*B*R)
				if (operB.getExists()){
					Utilities.subEquals(Temp,Utilities.mult(blockVectorY,
						Utilities.mult(Utilities.trans(blockVectorBY),Temp)));
					Utilities.setMatrix(blockVectorR,Temp,mskIdx);
				}
				// R=R-Y*(Y^T*R)
				else  {
					Utilities.subEquals(Temp,Utilities.mult(blockVectorY,
						Utilities.mult(Utilities.trans(blockVectorY),Temp)));
					Utilities.setMatrix(blockVectorR,Temp,mskIdx);
				}
			}

			// Apply preconditioner operT to active residuals
			// R(index)=operT*R(index)
			if (operT.getExists()){
				Utilities.setMatrix(blockVectorR,
					operT.operatorAction(Utilities.getMatrix(blockVectorR,mskIdx)),mskIdx);
				// project R onto subspace B-orthogonal to Y (or orthogonal to Y)
				if (blockVectorY.numColumns() != 0){
					Temp=Utilities.getMatrix(blockVectorR,mskIdx);
					// R=R-Y*(Y^T*B*R)
					if (operB.getExists()){
						Utilities.subEquals(Temp,Utilities.mult(blockVectorY,
							Utilities.mult(Utilities.trans(blockVectorBY),Temp)));
						Utilities.setMatrix(blockVectorR,Temp,mskIdx);
					}
					// R=R-Y*(Y^T*R)
					else  {
						Utilities.subEquals(Temp,Utilities.mult(blockVectorY,
							Utilities.mult(Utilities.trans(blockVectorY),Temp)));
						Utilities.setMatrix(blockVectorR,Temp,mskIdx);
					}
				}
			}
			
			// Make active (preconditioned) residuals orthogonal (B-orthogonal) to X
			Temp=Utilities.getMatrix(blockVectorR,mskIdx);
			if (operB.getExists()){
				// R(index)=R(index) - X((BX)'*R(index))
				Utilities.subEquals(Temp,Utilities.mult(blockVectorX,
					Utilities.mult(Utilities.trans(blockVectorBX),Temp)));
				Utilities.setMatrix(blockVectorR,Temp,mskIdx);
				}
			else {
				// R(index)=R(index) - X((X'*R(index))
				Utilities.subEquals(Temp,Utilities.mult(blockVectorX,
					Utilities.mult(Utilities.trans(blockVectorX),Temp)));
				Utilities.setMatrix(blockVectorR,Temp,mskIdx);
			}
			
			// B-orthonormalize R
			if (operB.getExists()){
				// BR(mskIdx)=B * R(mskIdx)
				Utilities.setMatrix(blockVectorBR,
					operB.operatorAction(Utilities.getMatrix(blockVectorR,mskIdx)),mskIdx);
				// gramRBR=R(mskIdx)^T * BR(mskIdx)
				gramRBR=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorR,mskIdx)),
					Utilities.getMatrix(blockVectorBR,mskIdx));
			}
			// gramRBR=R(mskIdx)^T * BR(mskIdx)
			else gramRBR=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorR,mskIdx)),
					Utilities.getMatrix(blockVectorR,mskIdx));
			
			Utilities.symmetrizeMatrix(gramRBR);
			Utilities.setMatrix(blockVectorR,Utilities.setOrthCholesky(gramRBR,
				Utilities.getMatrix(blockVectorR,mskIdx)),mskIdx);
			if (operB.getExists())  
				Utilities.setMatrix(blockVectorBR,
				Utilities.setOrthCholesky(gramRBR,Utilities.getMatrix(blockVectorBR,mskIdx)),mskIdx);
			Utilities.setMatrix(blockVectorAR,
				operA.operatorAction(Utilities.getMatrix(blockVectorR,mskIdx)),mskIdx);
			
			if (iterationNumber>1 && restartFlag==false){
				// B-orthonormalize P
				if (operB.getExists()) 
					gramPBP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorP,mskIdx)),
					Utilities.getMatrix(blockVectorBP,mskIdx));
				else gramPBP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorP,mskIdx)),
					Utilities.getMatrix(blockVectorP,mskIdx));
				
				Utilities.symmetrizeMatrix(gramPBP);
				Utilities.setMatrix(blockVectorP,Utilities.setOrthCholesky(gramPBP,
					Utilities.getMatrix(blockVectorP,mskIdx)),mskIdx);
				Utilities.setMatrix(blockVectorAP,Utilities.setOrthCholesky(gramPBP,
					Utilities.getMatrix(blockVectorAP,mskIdx)),mskIdx);
				if (operB.getExists()) Utilities.setMatrix(blockVectorBP,
					Utilities.setOrthCholesky(gramPBP,Utilities.getMatrix(blockVectorBP,mskIdx)),mskIdx);
						
				// Compute gram matrices
				gramXAR=Utilities.mult(Utilities.trans(blockVectorAX),Utilities.getMatrix(blockVectorR, mskIdx));
				gramXAP=Utilities.mult(Utilities.trans(blockVectorAX),Utilities.getMatrix(blockVectorP, mskIdx));
				gramRAR=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorAR, mskIdx)),
					Utilities.getMatrix(blockVectorR, mskIdx));
				Utilities.symmetrizeMatrix(gramRAR);
				
				gramRAP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorR, mskIdx)),
						Utilities.getMatrix(blockVectorAP, mskIdx));
				
				gramPAP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorP, mskIdx)),
						Utilities.getMatrix(blockVectorAP, mskIdx));
				Utilities.symmetrizeMatrix(gramPAP);
				
				if (operB.getExists()){
					gramXBR=Utilities.mult(Utilities.trans(blockVectorX),Utilities.getMatrix(blockVectorBR, mskIdx));
					gramXBP=Utilities.mult(Utilities.trans(blockVectorX),Utilities.getMatrix(blockVectorBP, mskIdx));
					gramRBP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorR, mskIdx)),
						Utilities.getMatrix(blockVectorBP, mskIdx));
					}
				else {
					gramXBR=Utilities.mult(Utilities.trans(blockVectorX),Utilities.getMatrix(blockVectorR, mskIdx));
					gramXBP=Utilities.mult(Utilities.trans(blockVectorX),Utilities.getMatrix(blockVectorP, mskIdx));
					gramRBP=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorR, mskIdx)),
						Utilities.getMatrix(blockVectorP, mskIdx));
				}
				
				// Build gramA matrix
				gramA=Utilities.blockMatrix(3,3,EIG,gramXAR,gramXAP,
					Utilities.trans(gramXAR),gramRAR,gramRAP,
					Utilities.trans(gramXAP),Utilities.trans(gramRAP),gramPAP);
					
				// Build gramB matrix
				nn=gramRBP.numRows();
				gramB=Utilities.blockMatrix(3,3,Matrices.identity(blockSize),null,gramXBP,
						null,Matrices.identity(nn),gramRBP,
						Utilities.trans(gramXBP),Utilities.trans(gramRBP),Matrices.identity(nn));
				
			} else {
				// Handle first iteration to compute gram matrices
				gramXAR=Utilities.mult(Utilities.trans(blockVectorAX),Utilities.getMatrix(blockVectorR, mskIdx));
				gramRAR=Utilities.mult(Utilities.trans(Utilities.getMatrix(blockVectorAR, mskIdx)),
					Utilities.getMatrix(blockVectorR, mskIdx));
				Utilities.symmetrizeMatrix(gramRAR);
				
				if (operB.getExists()) gramXBR=Utilities.mult(Utilities.trans(blockVectorX),
					Utilities.getMatrix(blockVectorBR, mskIdx));
				else gramXBR=Utilities.mult(Utilities.trans(blockVectorX),Utilities.getMatrix(blockVectorR, mskIdx));

				gramA=Utilities.blockMatrix(2,2,EIG,gramXAR,
					Utilities.trans(gramXAR),gramRAR);
				
				// Build upper triangular part of gramB matrix
				gramB=(DenseMatrix) Matrices.identity(blockSize+gramXBR.numColumns());
			}

			currentSubspaceDimension=gramA.numColumns();
			eigenvalues=new double[gramA.numColumns()];
			
			// Check condition number for restart
			conditionG=Math.log10(Utilities.getConditionNumber(gramB.copy()))+1;
			k=Math.max(0,iterationNumber -11 -(int)Math.round(Math.log(currentBlockSize)));
			conditionGMean=0;
			for (int i=k;i<iterationNumber;++i) conditionGMean=conditionGMean+conditionGHistory[i];
			conditionGMean=conditionGMean/(iterationNumber-k);
			
			if (eigensolverDebug) {
				conditionGramA=Utilities.getConditionNumber(gramA);
				conditionGramB=Utilities.getConditionNumber(gramB);
			}
			
			if (Utilities.eig(gramA,gramB,eigenvalues)!=0){
				System.out.println("Generalized eigenvalue solution could not be completed (LAPACK - DSYGV)");
			}
			
			// Assign top blockSize eigenvalues to diagonal of EIG
			EIG=new DenseMatrix(blockSize,blockSize);
			for (int i=0;i<blockSize;++i){
				EIG.set(i,i,eigenvalues[i]);
			}
			
			// Keep first blockSize columns of gramA
			if (iterationNumber>1){
				int i0=0;
				int i1=blockVectorX.numColumns()-1;
				int j0=0;
				int j1=blockVectorX.numColumns()-1;
				CX=Utilities.getMatrix(gramA,i0,i1,j0,j1);
				i0=i1+1;
				i1 += currentBlockSize;
				CR=Utilities.getMatrix(gramA,i0,i1,j0,j1);
				i0=i1+1;
				i1 += currentBlockSize;
				CP=Utilities.getMatrix(gramA,i0,i1,j0,j1);

				// Update matrices
				if (restartFlag){
					// Use block steepest descent to regain stability
					blockVectorP=Utilities.mult(Utilities.getMatrix(blockVectorR, mskIdx), CR);
					blockVectorAP=Utilities.mult(Utilities.getMatrix(blockVectorAR, mskIdx), CR);
					if (operB.getExists()) blockVectorBP=Utilities.mult(Utilities.getMatrix(blockVectorBR, mskIdx),CR);
				} else {
					blockVectorP=Utilities.add(Utilities.mult(Utilities.getMatrix(blockVectorR, mskIdx), CR),
						Utilities.mult(Utilities.getMatrix(blockVectorP, mskIdx), CP));
					blockVectorAP=Utilities.add(Utilities.mult(Utilities.getMatrix(blockVectorAR, mskIdx), CR),
						Utilities.mult(Utilities.getMatrix(blockVectorAP, mskIdx), CP));
					if (operB.getExists()) 
						blockVectorBP=Utilities.add(Utilities.mult(Utilities.getMatrix(blockVectorBR, mskIdx), CR),
						Utilities.mult(Utilities.getMatrix(blockVectorBP, mskIdx), CP));
				}
				
				blockVectorX=Utilities.add(Utilities.mult(blockVectorX, CX),blockVectorP);
				blockVectorAX=Utilities.add(Utilities.mult(blockVectorAX, CX),blockVectorAP);
				if (operB.getExists()) blockVectorBX=Utilities.add(Utilities.mult(blockVectorBX, CX),blockVectorBP);
			} else {
				// Obtain rows of gramA according to column widths of X and W
				int i0=0;
				int i1=blockVectorX.numColumns()-1;
				int j0=0;
				int j1=blockVectorX.numColumns()-1;
				CX=Utilities.getMatrix(gramA,i0,i1,j0,j1);
				i0=i1+1;
				i1 += blockVectorR.numColumns();
				CR=Utilities.getMatrix(gramA,i0,i1,j0,j1);
					
				blockVectorP=Utilities.mult(Utilities.getMatrix(blockVectorR, mskIdx), CR);
				blockVectorAP=Utilities.mult(Utilities.getMatrix(blockVectorAR, mskIdx), CR);
				blockVectorX=Utilities.add(Utilities.mult(blockVectorX, CX),
					Utilities.getMatrix(blockVectorP, mskIdx));
				blockVectorAX=Utilities.add(Utilities.mult(blockVectorAX,CX),blockVectorAP);
				if (operB.getExists()){
					blockVectorBP=Utilities.mult(Utilities.getMatrix(blockVectorBR,mskIdx),CR);
					blockVectorBX=Utilities.add(Utilities.mult(Utilities.getMatrix(blockVectorBX, mskIdx),CX),
						Utilities.getMatrix(blockVectorBP, mskIdx));
				}
			}
			
			// Compute residuals
			if (operB.getExists()) blockVectorR=Utilities.sub(blockVectorAX, Utilities.mult(blockVectorBX, EIG));
			else blockVectorR=Utilities.sub(blockVectorAX, Utilities.mult(blockVectorX, EIG));
			
			// Update active mask and active mask array
			for (int j=0;j<blockSize;++j){
				residualNorms[j]=Utilities.getMatrix(blockVectorR,j,j).norm(Norm.Frobenius);
				if (iterationNumber>1){ 
					if (residualNorms[j]<residualTolerance){
						// Remove index if in list
						if (activeMask.indexOf(j)>-1) activeMask.remove(activeMask.indexOf(j));
					}
				}
			}
			
			// Update active mask
			tempInteger=new Integer[activeMask.size()];
			tempInteger=activeMask.toArray(tempInteger);
			mskIdx=new int [activeMask.size()];
			for (int i=0; i < tempInteger.length; i++) mskIdx[i]=tempInteger[i].intValue();
			
			// Load data for history data collection
			for (int i=0;i<blockSize;++i){
				eigenvalueHistory.set(i, iterationNumber-1, eigenvalues[i]);
				residualNormsHistory.set(i, iterationNumber-1, residualNorms[i]);
			}
			residualMaximumHistory[iterationNumber-1]=getMaxResidualNorm();
			timeHistory[iterationCount-1]=System.currentTimeMillis()-timeIteration;
			if (Utilities.getMemoryUsage()>maxMemoryUsed) maxMemoryUsed=Utilities.getMemoryUsage();
			conditionGHistory[iterationCount]=conditionG;
			
			// Print out information concerning progress and convergence
			if (verbosityLevel>0){
				for (int i=0;i<blockSize;++i){
					System.out.printf("Eigenvalue lambda %17.16e \n",eigenvalues[i]);
				}
		        for (int i=0;i<blockSize;++i){
					System.out.printf("Residual norm %e \n",residualNorms[i]);
				}
		        System.out.printf("Iteration number: %d, Max iterations: %d\n",
	        		iterationNumber,maxIterations);
		        if (verbosityLevel>1){
		        	System.out.printf("Maximum residual:           \t%e (residual tol = %.2e)\n",
		        			getMaxResidualNorm(),residualTolerance);
		        	System.out.printf("Current block size:         \t%d \n",currentBlockSize);
		        	System.out.printf("Current subspace dimension: \t%d\n", currentSubspaceDimension);
		        	if (informationString.length()>0) System.out.println(informationString);
		        }
		        System.out.printf("\n");
			}
			
			if (verbosityLevel==-1)
	        	System.out.printf("Iteration number: %d Maximum residual: %e (residual tol = %.2e)\n",
	        	iterationNumber,getMaxResidualNorm(),residualTolerance);
			
			if (eigensolverDebug) 
				debug(conditionGramA,conditionGramB,blockVectorX,blockVectorR,blockVectorP,
				blockVectorBX,blockVectorBR,blockVectorBP,mskIdx);
		}
	
		// Save number of iterations that were executed
		setNumbIterationsExecuted(iterationCount);
		
		// Compact eigenvalue, residual norms and time data
		eigenvalueHistory=Utilities.getMatrix(eigenvalueHistory,0,blockSize-1,0,iterationCount-1);
		residualNormsHistory=Utilities.getMatrix(residualNormsHistory,0,blockSize-1,0,iterationCount-1);
		residualMaximumHistory[iterationCount-1]=getMaxResidualNorm();

		double[] tempCompact=new double[iterationCount];
		for (int i=0;i<iterationCount;++i) tempCompact[i]=residualMaximumHistory[i];
		residualMaximumHistory=tempCompact;
		
		// Save actual condition numbers
		tempCompact=new double[iterationCount+1];
		tempCompact[0]=Math.pow(conditionGHistory[0],10);
		for (int i=1;i<iterationCount+1;++i) tempCompact[i]=(Math.pow(conditionGHistory[i]-1,10));
		conditionGHistory=tempCompact;
		
		long[] tempTime=new long[iterationCount];
		for (int i=0;i<iterationCount;++i)tempTime[i]=timeHistory[i];
		timeHistory=tempTime;

		// Print out final information
		long t1 = System.currentTimeMillis();
		totalTime=t1-t0;
		if (verbosityLevel>0) {
			if (convergenceFlag) System.out.println("Lobpcg has converged");
			else System.out.println("Lobpcg has not converged to residual tolerance");
			System.out.printf("Number of iterations:       \t%d \n",iterationCount);
			System.out.printf("Final maximum residual norm:\t%e\n",getMaxResidualNorm());
			if (verbosityLevel>1) {
				System.out.printf("Execution time:             \t%-10.5f seconds (%d ms)\n",(t1-t0)/1000.,t1-t0);
				System.out.printf("Seconds/iteration:          \t%-10.5f seconds (%d ms)\n",
					(t1-t0)/(1000.*iterationCount),(t1-t0)/iterationCount);
				System.out.printf("Maximum memory used:        \t%.2fM\n",maxMemoryUsed/1000000.);
			}
			if (restartCount>0) System.out.printf("Number of restarts:        \t%d\n",restartCount);
			System.out.println("Exiting Lobpcg!\n\n");
		}
		if (convergenceFlag) return(1); // Converged
		else return(-1);  // Did not converge
	}
	
	// Print out debug information
	void debug(double conditionGramA,double conditionGramB,DenseMatrix blockVectorX,
		DenseMatrix blockVectorR,DenseMatrix blockVectorP,DenseMatrix blockVectorBX,
		DenseMatrix blockVectorBR,DenseMatrix blockVectorBP,int[] mskIdx){
		DenseMatrix Temp;
		DenseMatrix Temp1;
		DenseMatrix Temp2;
		
		// Check condition numbers of gramA and gramB
		System.out.printf("conditionGramA=%e\n",conditionGramA);
		System.out.printf("conditionGramB=%e\n",conditionGramB);
		// Check ranks
		System.out.printf("rank blockVectorX = %d \n",Utilities.getRank(blockVectorX));
		System.out.printf("rank blockVectorR = %d \n",Utilities.getRank(blockVectorR));
		System.out.printf("rank blockVectorP = %d \n",Utilities.getRank(blockVectorP));

		if (operB.getExists()){
			System.out.printf("norm(X'*B*X - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(blockVectorX),
				blockVectorBX),Matrices.identity(blockVectorX.numColumns())).norm(Norm.Frobenius));
			Temp1=Utilities.getMatrix(blockVectorR,mskIdx);
			Temp2=Utilities.getMatrix(blockVectorBR,mskIdx);
			System.out.printf("norm(R(idx)'*B*R(idx) - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(Temp1),
				Temp2),Matrices.identity(Temp1.numColumns())).norm(Norm.Frobenius));
			Temp1=Utilities.getMatrix(blockVectorP,mskIdx);
			Temp2=Utilities.getMatrix(blockVectorBP,mskIdx);
			System.out.printf("norm(P(idx)'*B*P(idx) - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(Temp1),
				Temp2),Matrices.identity(Temp1.numColumns())).norm(Norm.Frobenius));
		}
		else {
			System.out.printf("norm(X'*X - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(blockVectorX),
				blockVectorX),Matrices.identity(blockVectorX.numColumns()))
					.norm(Norm.Frobenius));
			Temp=Utilities.getMatrix(blockVectorR,mskIdx);
			System.out.printf("norm(R(idx)'*R(idx) - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(Temp),Temp),
				Matrices.identity(Temp.numColumns())).norm(Norm.Frobenius));
			System.out.printf("norm(R'*R - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(blockVectorR),blockVectorR),
				Matrices.identity(blockVectorR.numColumns())).norm(Norm.Frobenius));
			System.out.printf("norm(X'*R)=%e\n", 
				Utilities.mult(Utilities.trans(blockVectorX),blockVectorR).norm(Norm.Frobenius));
			Temp=Utilities.getMatrix(blockVectorP,mskIdx);
			System.out.printf("norm(P(idx)'*P(idx) - I)=%e\n", 
				Utilities.sub(Utilities.mult(Utilities.trans(Temp),
				Temp),Matrices.identity(Temp.numColumns()))
					.norm(Norm.Frobenius));
		}
		System.out.printf("Norm of residuals = %e\n",blockVectorR.norm(Norm.Frobenius));
		// Check to see if X is orthogonal to Y
		if (blockVectorY.numColumns() != 0){
			if (operB.getExists()) 
				System.out.printf("Checking if constraints are orthogonal to X norm(Y^T*B*X)=%e\n",
					Utilities.mult(Utilities.trans(blockVectorY),blockVectorBX).norm(Norm.Frobenius));
			else System.out.printf("Checking if constraints are orthogonal to X norm(Y^T*X)=%e\n",
					Utilities.mult(Utilities.trans(blockVectorY),blockVectorX).norm(Norm.Frobenius));
			Utilities.printMatrixInformation(blockVectorX);
			Utilities.printMatrixInformation(blockVectorY);
		}
		// Print time for last matrix multiply
		System.out.printf("Time for operA: %10.5f seconds (%d ms)\n",
			operA.getTimeOperatorApply()/1000.,operA.getTimeOperatorApply());
		System.out.printf("Time for operB: %10.5f seconds (%d ms)\n",
				operB.getTimeOperatorApply()/1000.,operB.getTimeOperatorApply());
		System.out.printf("Time for operT: %10.5f seconds (%d ms)\n",
				operT.getTimeOperatorApply()/1000.,operT.getTimeOperatorApply());
		System.out.printf("\n");
	}
}
