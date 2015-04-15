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
import no.uib.cipr.matrix.DenseCholesky;
import no.uib.cipr.matrix.LowerSPDDenseMatrix;

/**
SEJ Java OperatorPrecDenseCholesky class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
@author Rico Argentati
*/
public class OperatorPrecDenseCholesky extends Operator {
	
	DenseCholesky ch;
	/* ------------------------
	   Constructors
	 * ------------------------ */
	public OperatorPrecDenseCholesky(){super();} // Must have a default constructor
	public OperatorPrecDenseCholesky(DenseMatrix T){ 
		ch = new DenseCholesky(T.numColumns(),false);
		LowerSPDDenseMatrix LL= new LowerSPDDenseMatrix(T);
		ch.factor(LL);
		setOperatorSize(T.numColumns());
	}

	// This overrides operatorAction in Operator class
	public DenseMatrix operatorAction(DenseMatrix X){
		long t = System.currentTimeMillis();
		// check compatibility of matrix dimensions
		if (getOperatorSize()!=X.numRows()){
			System.out.println("Preconditioner not compatible with input matrix size, "+
					getOperatorSize()+" "+X.numRows());
			return (null);
		}
		ch.solve(X);  // X=inv(T)*X
		timeOperatorApply=System.currentTimeMillis()-t;
		return X;
	}
}
