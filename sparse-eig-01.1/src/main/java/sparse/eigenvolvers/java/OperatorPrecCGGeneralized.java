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
	import no.uib.cipr.matrix.DenseVector;
	import no.uib.cipr.matrix.MatrixEntry;
	import no.uib.cipr.matrix.Vector;
	import no.uib.cipr.matrix.sparse.AbstractIterationMonitor;
	import no.uib.cipr.matrix.sparse.CG;
	import no.uib.cipr.matrix.sparse.CompRowMatrix;
	import no.uib.cipr.matrix.sparse.IterativeSolver;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;

	/**
	SEJ Java public class OperatorPrecCGGeneralized class.
	<p>
	This class is experimental
	<p>
	This class sets up a conjugate gradient solver using a shift. The goal is to multiply the
	input block matrix, using {@link #operatorAction(DenseMatrix)}, 
	by an approximate inverse to the operator T
	<p>
	SEJ uses the MTJ Java library (matrix-toolkits-java)
	and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
	@author Rico Argentati
	*/
	public class OperatorPrecCGGeneralized extends EigPreconditioner {
		/* ------------------------
		   Class variables
		 * ------------------------ */
		CompRowMatrix A;
		CompRowMatrix B;
			
		/* ------------------------
		   Constructors
		 * ------------------------ */
		public OperatorPrecCGGeneralized(){super();} // Must have a default constructor
		
		/**
	     * Constructor sets up conjugate gradient solver using operator T 
	     * in compressed row format
	     * 
	     * @param T compressed row matrix
	     *           
	     */
		public OperatorPrecCGGeneralized(CompRowMatrix A,CompRowMatrix B,SparseEigensolver eigSolv){
			this.A=A.copy();
			this.B=B.copy();
			setOperatorSize(A.numColumns());
			this.eigSolv=eigSolv;  // need calling solver to get smallest active Ritz value
			operatorSparse=true;
			currentShift=1;
			
			// Setup solver
			DenseVector v=new DenseVector(getOperatorSize());
			solv=new CG(v);
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
			
			// Construct precondition matrix operatorSparseMatrix = A - \mu * B
			// where \mu is the smallest active residual
			double shiftRatio=eigSolv.getSmallestActiveRitzValue()/currentShift;
			for (MatrixEntry e : B)  e.set(e.get()*shiftRatio); 
			operatorSparseMatrix=Utilities.sub(A, B);
			currentShift=eigSolv.getSmallestActiveRitzValue();
			
			// Solve using CG
			solv.setIterationMonitor(new SimpleIterationMonitor(cgIterations)); // set number of iterations
			DenseVector v=new DenseVector(X.numRows());
			for (int j=0; j<X.numColumns(); ++j){
				v=Utilities.getVectorMatrix(X,j);
				try {
					solv.solve(operatorSparseMatrix, v, v);
				}
				catch (IterativeSolverNotConvergedException e){}
				Utilities.setVectorMatrix(X,v,j);
			}
			return X;
		}

		// This class is need to input maximum iterations to solver
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
