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

/**
SEJ Java SparseEigensolverConstants class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
@author Rico Argentati
*/
public class SparseEigensolverConstants {
	public static final int CG_NUMBER_ITERATIONS=10;    // used for eigenvalue precondioning solver 
	public static final double RESIDUAL_TOLERANCE=1e-5; // not used since actually computed
	public static final int MAX_NUMBER_ITERATIONS=20;	// default maximum number of iterations
	public static final int VERBOSITY_LEVEL=1;          // print out basic information at each iteration
	public static final int BLOCK_SIZE=10;				// default block size
	public static final int RESTART_COUNT_BLOCKDAVIDSON=5;
	public static final String PRINT_STRING="%17.15e ";	// defaults print string for print utilities
	// TMJ solver preconditioners
	public enum subPRECONDIONERS {DIAGONAL, SSOR, ICC, ILU, ILUT, AMG}
}
