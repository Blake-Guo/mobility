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

import no.uib.cipr.matrix.sparse.*;

/**
SEJ Java EigPreconditioner class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
This class provides the setup for the solver and TMJ preconditioning

@author Rico Argentati
*/

public class EigPreconditioner extends Operator {
	/* ------------------------
	   Class variables
	 * ------------------------ */
	CompRowMatrix T; // Matrix input for solver
	int cgIterations=SparseEigensolverConstants.CG_NUMBER_ITERATIONS; // Inner iteration default
	IterativeSolver solv; 		// Iterative TMJ solver
	SparseEigensolver eigSolv;  // Sparse eigensolver (e.g. Lobpcg)
	double currentShift;		// Used for shift in solver

	/**
     * Check input for TMJ preconditiner for solver 
     * 
     * @param s name of TMJ preconditioner (DIAGONAL, SSOR, ICC, ILU, ILUT, AMG)
     *           
     */
	public boolean checkSubPreconditionerName(String s){
		for (SparseEigensolverConstants.subPRECONDIONERS prec :
			SparseEigensolverConstants.subPRECONDIONERS.values()) { 
			if (prec.name().equals(s)){
				return true;
			}
		} 
		return false;
	}
	
	/**
     * Activate sub TMJ sub-preconditioner for solver 
     * 
     * @param T compressed row matrix (inv(T) approximates matrix used to solve eigenvalue problem
     * @param solv iterative solver
     * @param prec name of sub-preconditioner (DIAGONAL, SSOR, ICC, ILU, ILUT, AMG)
     * 
     * @return 1 sub-preconditioner assigned, -1 sub-preconditioner not assigned,
     *           
     */
	 protected int activateSubPreconditioner(CompRowMatrix T,IterativeSolver solv, String prec){
		Preconditioner M=solv.getPreconditioner();
	    if (prec.equals(SparseEigensolverConstants.subPRECONDIONERS.DIAGONAL.toString())){
	    	M = new DiagonalPreconditioner(T.numRows());
	    	System.out.println("Got DIAGONAL TMJ sub-preconditioner");
	    }
	    else if (prec.equals(SparseEigensolverConstants.subPRECONDIONERS.SSOR.toString())){
	    	M = new SSOR(new CompRowMatrix(T));
	    	System.out.println("Got SSOR TMJ sub-preconditioner");
	    }
	    else if (prec.equals(SparseEigensolverConstants.subPRECONDIONERS.ICC.toString())){
	    	M = new ICC(new CompRowMatrix(T));
	    	System.out.println("Got ICC TMJ sub-preconditioner");
	    }
	    else if (prec.equals(SparseEigensolverConstants.subPRECONDIONERS.ILU.toString())){
	    	M = new ILU(new CompRowMatrix(T));
	    	System.out.println("Got ILU TMJ sub-preconditioner");
	    }  
	    else if (prec.equals(SparseEigensolverConstants.subPRECONDIONERS.AMG.toString())){
	    	M = new AMG();
	    	System.out.println("Got AMG TMJ sub-preconditioner");
	    }  
	    else return -1;
	    M.setMatrix(T);
	    // Attach the sub-preconditioner
	    solv.setPreconditioner(M);
	    return 1;
	}
	
}
