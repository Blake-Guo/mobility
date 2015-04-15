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

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import no.uib.cipr.matrix.sparse.CompRowMatrix; 

/**
SEJ Java Operator class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
The Operator class provides convenient methods for defining an operator, in a flexible 
manner,  that implements a matrix block vector multiply.
A user can implement whatever operation is desired desire by
extending this class and overriding {@link #operatorAction}.

@author Rico Argentati
*/
public class Operator {
	/* ------------------------
	   Class variables
	 * ------------------------ */
	protected DenseMatrix operatorMatrix;
	protected CompRowMatrix operatorSparseMatrix; 
	protected boolean operatorSparse=false;
	private int operatorSize=0;
	protected long timeOperatorApply=0;
	
	/* ------------------------
	   Constructors
	 * ------------------------ */
	public Operator(){}
	
	/**
     * Constructor sets up operator using dense matrix A
     * 
     * @param A dense matrix
     *           
     */
	public Operator(DenseMatrix A){
		operatorMatrix=A;
		setOperatorSize(A.numColumns());
		operatorSparse=false;
	}
	
	/**
     * Constructor sets up operator using a compressed row matrix A
     * 
     * @param A compressed row matrix
     *           
     */
	public Operator(CompRowMatrix A){
		operatorSparseMatrix=A;
		setOperatorSize(A.numColumns());
		operatorSparse=true;
	}
	
	/* ------------------------
	   Public Methods
	 * ------------------------ */
	/**
     * get operator size (dimension)
     * 
     * @return operator size (dimensions)
     *           
     */
	public int getOperatorSize(){
		return operatorSize;
	}
	
	/**
     * Check if operator exists
     * 
     * @return <code>true</code> or <code>false</code>
     *           
     */
	public boolean getExists(){
		if (operatorSize>0) return true;
		else return false;
	}
	
	/**
     * Set operator to identity
     * 
     * @param n size of identity (n x n)
     *           
     */
	public void setOperatorIdentity(int n){
		operatorMatrix=(DenseMatrix) Matrices.identity(n);
		setOperatorSize(operatorMatrix.numColumns());
		operatorSparse=false;
	}
	
	/**
     * Set operator to dense matrix
     * 
     * @param A dense matrix
     *           
     */
	public void setOperatorLoadMatrix(DenseMatrix A){
		operatorMatrix=A;
		setOperatorSize(operatorMatrix.numColumns());
		operatorSparse=false;
	}
	
	/**
     * Set operator to compressed row matrix
     * 
     * @param A compressed row matrix
     *           
     */
	public void setOperatorLoadMatrix(CompRowMatrix A){
		operatorSparseMatrix=A;
		setOperatorSize(A.numColumns());
		operatorSparse=true;
	}
	
	/**
     * Set operator to dense random symmetric matrix
     * 
     * @param n size of matrix (n x n)
     *           
     */
	public void setOperatorRandomSymmetric(int n){
		operatorMatrix=(DenseMatrix) Matrices.random(n,n);
		//Make it symmetric
		Utilities.addEquals(operatorMatrix, Utilities.trans(operatorMatrix));
		setOperatorSize(operatorMatrix.numColumns());
		operatorSparse=false;
	}
	
	/**
     * Set operator to a dense random positive definite matrix
     * 
     * @param n size of matrix (n x n)
     *           
     */
	public void setOperatorRandomPositiveDef(int n){
		operatorMatrix=(DenseMatrix) Matrices.random(n,n);
		//Make it symmetric and hopefully positive definite 
		operatorMatrix=Utilities.mult(operatorMatrix, Utilities.trans(operatorMatrix));
		setOperatorSize(operatorMatrix.numColumns());
		operatorSparse=false;
	}
	
	/**
     * Get dense matrix representing operator. Caution: operator must be small to not run out of memory. 
     * This is primarily for testing
     * 
     * @return dense matrix
     *           
     */
	public DenseMatrix getDenseOperatorMatrix(){
		return (operatorAction(Matrices.identity(getOperatorSize())));
	}
	
	/**
     * Get time for last operator action/multiply
     * 
     * @return time in ms
     *           
     */
	public long getTimeOperatorApply(){
		return timeOperatorApply;
	}
	
	/**
     * This method applies the operator to a dense input matrix
     * and returns a dense matrix. The method {@link #operatorAction} can be overwritten
     * to provide any matrix block vector multiply or any operation desired by the user
     * as long as <code>X</code> is returned.
     * 
     * @param X matrix
     * 
     * @return <code>X</code> matrix
     *           
     */
	public DenseMatrix operatorAction(DenseMatrix X){
		long t = System.currentTimeMillis();
		// check compatibility of matrix dimensions
		if (getOperatorSize()!=X.numRows()){
			System.out.println("Operator not compatible with input matrix size, "+
					getOperatorSize()+" "+X.numColumns());
			return (null);
		}
		if (operatorSparse){
			DenseMatrix Temp=new DenseMatrix(X.numRows(),X.numColumns());
			operatorSparseMatrix.mult(X,Temp); //X=operatorSparseMatrix*X
			timeOperatorApply=System.currentTimeMillis()-t;
			return(Temp);
		} else {
			X=Utilities.mult(operatorMatrix, X);
			timeOperatorApply=System.currentTimeMillis()-t;
			return X;
		}
	}
	
	/* ------------------------
	   Protected Methods
	 * ------------------------ */
	
	/**
     * Set operator size (dimension)
     * 
     * @param operatorSize (operatorSize x operatorSize)
     *           
     */
	protected void setOperatorSize(int operatorSize){
		this.operatorSize=operatorSize;
	}
}
