/*
 * Copyright (C) 2012 Rico Argentati
 * 
 * This file is part of SEJ (Sparse Eigensolvers for Java).
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General protected License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General protected License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General protected License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package sparse.eigenvolvers.java;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.sparse.CompRowMatrix;

/**
SEJ Java SparseEigensolver class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
The SparseEigensolver abstract class provides class variables and methods that can be used
to implement a general sparse eigenvalue solver.

@author Rico Argentati
*/
public abstract class SparseEigensolver {
 	
	/* ------------------------
	   Class variables**
	 * ------------------------ */
	protected DenseMatrix blockVectorX=new DenseMatrix(0,0); 	// Input block vector
	protected DenseMatrix blockVectorY=new DenseMatrix(0,0);	// Matrix of constraints
	protected int blockSize=SparseEigensolverConstants.BLOCK_SIZE;		// block size (column dimension of blockVectorX)
	protected int currentBlockSize;					// Current block size
	protected int maxIterations;					// Maximum iterations
	protected int numbIterationsExecuted;			// Number of iteration executed
	protected boolean maxIterationsSet=false;
	protected double residualTolerance;				// Residual tolerance
	protected boolean residualToleranceSet=false;
	protected String informationString=""; // Information string to display
	
	// Verbosity level: 0 no output, 1 detailed output,
	// >1 detailed output but only display on every verbosityLevel iteration count
	protected int verbosityLevel=SparseEigensolverConstants.VERBOSITY_LEVEL; 	
	protected double[] eigenvalues;   			// Vector of eigenvalues
	protected double[] residualNorms;			// Vector of residual norms
	protected long[] timeHistory;				// Store time in milliseconds for each iteration
	protected long totalTime=0;					// Total execution time in milliseconds
	protected DenseMatrix eigenvalueHistory;	// Eigenvalue history (blockSize x maxIterations)
	protected DenseMatrix residualNormsHistory;	// Residual norm history (blockSize x maxIterations)
	protected double[] residualMaximumHistory;  // Residual maximum history
	protected double[] conditionGHistory;  		// Condition number history (see lobpcg.m Matlab)
	protected Operator operA=new Operator();	// Main matrix operator 
	protected Operator operB=new Operator();	// Additional operator for generalized eigenvalue problem
	protected Operator operT=new Operator();	// Precondioner operator
	protected boolean eigensolverDebug=false;
	DenseMatrix saveBlockVectorX=null;			// Used to store location of blockVectorX so it can be returned
	
	/* ------------------------
	   protected Methods
	 * ------------------------ */
	/* ------------------------
	   Set and get methods
	 * ------------------------ */
	
	 /**
     * Set up primary operator 
     * 
     * @param operA
     *            Operator
     */
	public void setA(Operator operA){
		this.operA=operA;
	}
	
	 /**
     * Set up primary operator
     * 
     * @param A
     *            Matrix
     */
	public void setA(DenseMatrix A){
		this.operA=new Operator(A);
	}
	
	 /**
     * Set up primary operator
     * 
     * @param A
     *            CompRowMatrix
     */
	public void setA(CompRowMatrix A){
		this.operA=new Operator(A);
	}
	
	 /**
     * Set up secondary operator
     * 
     * @param operB
     *            Operator
     */
	public void setB(Operator operB){
		this.operB=operB;
	}
	
	/**
     * Set up secondary operator
     * 
     * @param B
     *            Matrix
     */
	public void setB(DenseMatrix B){
		this.operB=new Operator(B);
	}
	
	/**
     * Set up secondary operator
     * 
     * @param B
     *            CompRowMatrix
     */
	public void setB(CompRowMatrix B){
		this.operB=new Operator(B);
	}
	
	 /**
     * Set up preconditioner
     * 
     * @param operT
     *            Operator
     */
	public void setT(Operator operT){
		this.operT=operT;
	}
	
	 /**
     * Set up preconditioner
     * 
     * @param T
     *            Matrix
     */
	public void setT(DenseMatrix T){
		this.operT=new Operator(T);
	}
	
	 /**
     * Set up preconditioner
     * 
     * @param T
     *            CompRowMatrix
     */
	public void setT(CompRowMatrix T){
		this.operT=new Operator(T);
	}
	
	/**
     * Set up block vector to define initial subspace
     * 
     * @param blockVectorX
     *            Matrix
     */
	public void setBlockVectorX(DenseMatrix blockVectorX){
		this.blockVectorX=blockVectorX;
		if (this.blockVectorX.numColumns()>0) setBlockSize(this.blockVectorX.numColumns());
	}
	
	/**
     * Set block size of blockVectorX
     * <p>
     * Used when blockVectorX is set to random matrix
     * 
     * @param blockSize column dimension of blockVectorX
     */
	public void setBlockSize(int blockSize){
		this.blockSize=blockSize;
	}
	
	/**
     * Get blockSize 
     */
	public int getBlockSize(){
		return this.blockSize;
	}
	
	
	/**
     * Set numbIterationsExecuted
     * 
     * @param numbIterationsExecuted
     */
	protected void setNumbIterationsExecuted(int numbIterationsExecuted){
		this.numbIterationsExecuted=numbIterationsExecuted;
	}
	
	/**
     * Get numbIterationsExecuted
     */
	public int getNumbIterationsExecuted(){
		return(numbIterationsExecuted);
	}
	
	/**
     * Set up block vector to define contraints
     * 
     * @param blockVectorY
     *            matrix
     */
	public void setBlockVectorY(DenseMatrix blockVectorY){
		this.blockVectorY=blockVectorY;
	}

	/**
     * Set up maximum number of iterations
     * 
     * @param maxIterations
     *            
     */
	public void setMaxIterations(int maxIterations){
		this.maxIterations=maxIterations;
		maxIterationsSet=true;
	}
	
	/**
     * Set up residual tolerance
     * 
     * @param residualTolerance
     *            
     */
	public void setResidualTolerance(double residualTolerance){
		this.residualTolerance=residualTolerance;
		residualToleranceSet=true;
	}
	
	/**
     * Set up verbosity level
     * 
     * @param verbosityLevel
     *            =0 (no output)<br>
     *            =n>0 print information on progress every nth iteration<br>
     *            default = 1
     */
	public void setVerbosityLevel(int verbosityLevel){
		this.verbosityLevel=verbosityLevel;
	}
	
	/**
     * Get eigenvalues
     * 
     * @return 1D double array of eigenvalues
     *           
     */
	public double[] getEigenvalues(){
		double[] Temp=new double[blockVectorX.numColumns()];
		for (int i=0;i<blockVectorX.numColumns();++i) Temp[i]=eigenvalues[i];
		return(Temp);
	}
	
	/**
     * Get column matrix of eigenvalues
     * 
     * @return blockSize x 1 matrix of eigenvalues
     *           
     */
	public DenseMatrix getEigenvaluesMatrix(){
		DenseMatrix Z=new DenseMatrix(blockVectorX.numColumns(),1);
		for (int i=0;i<blockVectorX.numColumns();++i) Z.set(i,0,eigenvalues[i]);
		return Z;
	}
	
	/**
     * Get residuals norms
     * 
     * @return 1D double array of residual norms
     *           
     */
	public double[] getResidualNorms(){
		double[] Temp=new double[residualNorms.length];
		for (int i=0;i<residualNorms.length;++i) Temp[i]=residualNorms[i];
		return Temp;
	}
	
	/**
     * Get column matrix of residual norms
     * 
     * @return blockSize x 1 matrix of residual norms
     *           
     */
	public DenseMatrix getResidualNormsMatrix(){
		int nn=residualNorms.length;
		DenseMatrix Z=new DenseMatrix(nn,1);
		for (int i=0;i<nn;++i) Z.set(i,0,residualNorms[i]);
		return Z;
	}
	
	/**
     * Get block matrix of eigenvectors
     * 
     * @return n x blocksSize matrix of eigenvectors
     *           
     */
	public DenseMatrix getEigenvectors(){
		return(blockVectorX.copy());
	}
	
	/**
     * Get maximum residual norm
     * 
     * @return maximum residual norm
     *           
     */
	public double getMaxResidualNorm(){
	    double max = residualNorms[0];
	    for (int i = 1; i < residualNorms.length; i++) {
	          if (residualNorms[i] > max) max = residualNorms[i];
	    }
	    return max;
	}
	
	/**
     * Get eigenvalueHistory
     * 
     * @return 2D double array of eigenvalueHistory
     *           
     */
	public double[][] getEigenvalueHistory(){
		double[][] Temp=new double[eigenvalueHistory.numRows()][eigenvalueHistory.numColumns()];
		for (int i=0;i<eigenvalueHistory.numRows();i++) {
			 for (int j=0;j<eigenvalueHistory.numColumns();j++) {
				 Temp[i][j]=eigenvalueHistory.get(i,j);
			 }
		}
		return Temp;
	}
	
	/**
     * Get matrix of eigenvalueHistory
     * 
     * @return number of iterations x blockSize matrix of eigenvalueHistory
     *           
     */
	public DenseMatrix getEigenvalueHistoryMatrix(){
		return(eigenvalueHistory.copy());
	}
	
	
	/**
     * Get residualNormsHistory
     * 
     * @return 2D double array of residualNormsHistory
     *           
     */
	public double[][] getResidualNormsHistory(){
		double[][] Temp=new double[residualNormsHistory.numRows()][residualNormsHistory.numColumns()];
		for (int i=0;i<residualNormsHistory.numRows();i++) {
			 for (int j=0;j<residualNormsHistory.numColumns();j++) {
				 Temp[i][j]=residualNormsHistory.get(i,j);
			 }
		}
		return Temp;
	}
	
	/**
     * Get matrix of residualNormsHistory
     * 
     * @return number of iterations x blockSize matrix of residualNormsHistory
     *           
     */
	public DenseMatrix getResidualNormsHistoryMatrix(){
		return(residualNormsHistory.copy());
	}
	
	/**
     * Get residual norm maximum history
     * 
     * @return array of maximum norm hostory
     *           
     */
	public double[] getResidualMaximumHistory(){
		double[] Temp=new double[residualMaximumHistory.length];
		for (int i=0;i<residualMaximumHistory.length;++i) Temp[i]=residualMaximumHistory[i];
		return Temp;
	}
	
	/**
     * Get condition number history
     * 
     * @return array of condition number history
     *           
     */
	public double[] getConditionNumberHistory(){
		double[] Temp=new double[conditionGHistory.length];
		for (int i=0;i<conditionGHistory.length;++i) Temp[i]=conditionGHistory[i];
		return Temp;
	}
	
	/**
     * Get time history by iteration
     * 
     * @return array of time history (ms)
     *           
     */
	public long[] getTimeHistory(){
		long[] Temp=new long[timeHistory.length];
		for (int i=0;i<timeHistory.length;++i) Temp[i]=timeHistory[i];
		return Temp;
	}
	
	/**
     * Get total execution time
     * 
     * @return total execution time (ms)
     *           
     */
	public long getTotalTime(){
		return totalTime;
	}
	
	/**
     * Get the smallest active Ritz value. This is available for precondioning while an
     * eigensolver is running
     * 
     * @return smallest active Ritz value
     *           
     */
	public double getSmallestActiveRitzValue(){
		return eigenvalues[blockSize-currentBlockSize];
	}
	
	/**
     * Set print string to display for each iteration
     * 
     * @param s string
     *            
     */
	public void setInformationString(String s){
		informationString=s;
	}
	
	/**
     * Turn on debugging information output
     */
	public void setDebugOn(){
		eigensolverDebug=true;
	}

	/**
     * Turn off debugging information output (default)
     */
	public void setDebugOff(){
		eigensolverDebug=false;
	}
	
}
