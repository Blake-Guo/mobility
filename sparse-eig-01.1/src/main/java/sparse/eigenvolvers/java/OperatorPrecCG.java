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

//import sparse.eigenvolvers.java.EigPreconditioner.subPrecondioner;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.*;

/**
SEJ Java OperatorPrecCG class.
<p>
This class sets up a conjugate gradient solver. The goal is to multiply the
input block matrix, using {@link #operatorAction(DenseMatrix)}, 
by an approximate inverse to the operator T
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
@author Rico Argentati
*/
public class OperatorPrecCG extends EigPreconditioner {
	/* ------------------------
	   Constructors
	 * ------------------------ */
	public OperatorPrecCG(){super();} // Must have a default constructor
	
	/**
     * Constructor sets up conjugate gradient solver using operator T 
     * in compressed row format
     * 
     * @param T compressed row matrix
     *           
     */
	public OperatorPrecCG(CompRowMatrix T){
		this.T=T.copy();
		setOperatorSize(T.numColumns());
		
		// Setup solver
		DenseVector v=new DenseVector(getOperatorSize());
		solv=new CG(v);
	}
	
	public OperatorPrecCG(CompRowMatrix T, String prec){
		this.T=T.copy();
		setOperatorSize(T.numColumns());
		
		// Setup solver
		DenseVector v=new DenseVector(getOperatorSize());
		solv=new CG(v);
		
		// Setup sub-preconditioner
		if (checkSubPreconditionerName(prec)){
			activateSubPreconditioner(T,solv,prec);
		} else {
			System.out.println(getClass().getSimpleName()+" - invalid TMJ sub-preconditioner: "+prec);
			System.out.print("Use: ");
			for (SparseEigensolverConstants.subPRECONDIONERS p : 
				SparseEigensolverConstants.subPRECONDIONERS.values()) System.out.print((p.name()+" "));
			System.out.printf("\n");
			System.exit(0);
		}
	}
	
	/* ------------------------
	   Public Methods
	 * ------------------------ */
	/**
     * Set number of conjugate gradient (inner) iterations
     * 
     * @param cgIterations number of conjugate gradient (inner) iterations
     *           
     */
	public void setCGNumberIterations(int cgIterations){
		this.cgIterations=cgIterations;
	}
	
	/**
     * This method applies the approximate inverse of T to a dense input matrix
     * and returns a dense matrix. 
     * <p>
     * This method overrides {@link #operatorAction(DenseMatrix)}
     * 
     * @param X matrix
     * 
     * @return <code>X<code> matrix
     *           
     */
	public DenseMatrix operatorAction(DenseMatrix X){
		// check compatibility of matrix dimensions
		if (getOperatorSize()!=X.numRows()){
			System.out.println("Preconditioner not compatible with input matrix size, "+
					getOperatorSize()+" "+X.numRows());
			return (null);
		}
		// Solve using CG
		solv.setIterationMonitor(new SimpleIterationMonitor(cgIterations)); // set number of iterations
		DenseVector v=new DenseVector(X.numRows());
		for (int j=0; j<X.numColumns(); ++j){
			v=Utilities.getVectorMatrix(X,j);
			try {
				solv.solve(T, v, v);
			}
			catch (IterativeSolverNotConvergedException e){}
			Utilities.setVectorMatrix(X,v,j);
		}
		return X;
	}

	// This class is needed to input maximum iterations to solver
	private static class SimpleIterationMonitor extends AbstractIterationMonitor {
	    private int max;
	 
	    SimpleIterationMonitor(int max) {
	       this.max = max;
	     }
	    protected boolean convergedI(double r, Vector x) throws IterativeSolverNotConvergedException {
	       return convergedI(r);
	     }
	    protected boolean convergedI(double r) throws IterativeSolverNotConvergedException {
	       return iter >= max;
	     }
	}
}
