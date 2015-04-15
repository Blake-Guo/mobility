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
import no.uib.cipr.matrix.NotConvergedException;
import no.uib.cipr.matrix.SVD;

/**
SEJ Java Subspace class.
<p>
SEJ uses the MTJ Java library (matrix-toolkits-java)
and Netlib Java (netlib-java) for numerical linear algebra and matrix computations.
<p>
The class provides several methods for calculating angles between subspaces,
also known as principal angles. As one subspace becomes close to another
subspace the angles approach zero. 
There are min(dim(X),dim(Y)) principal angles. If all the angles are zero, then
either span{X} is a subspace of span{Y} or span{Y} is a subspace of span{X}.
Of course if X and Y have the same dimension and the angles are all zero, 
then the subspaces represented by the columns of X and Y are equal.
<br>
See <a href=http://en.wikipedia.org/wiki/Principal_angles<a>
<br>
See A. V. Knyazev and M. E. Argentati, Principal Angles between Subspaces 
in an A-Based Scalar Product: Algorithms and Perturbation Estimates. 
SIAM Journal on Scientific Computing, 23 (2002), no. 6, 2009-2041.
<a href=http://epubs.siam.org/sam-bin/dbq/article/37733<a>
<br>
Also see MATLAB program subspacea.m
<a href=http://www.mathworks.com/matlabcentral/fileexchange/55-subspacea-m<a>

@author Rico Argentati
*/

public class Subspace {
	
	// This method computes all principal angles using a sine based algorithm.
	// Small angles are computed more accurately than larger angles
	public static double[] getPrincipleAngles(DenseMatrix X,DenseMatrix Y){
		double[] s=new double[0];
		if (X.numColumns()==0 || Y.numColumns()==0 ) return s;
		DenseMatrix QX;
		DenseMatrix QY;
		
		// Orthonormalize X and Y
		if (X.numColumns()<Y.numColumns()){ //Switch X and Y in this case
			QX=Utilities.orth(Y);
			QY=Utilities.orth(X);
		} else {
			QX=Utilities.orth(X);
			QY=Utilities.orth(Y);
		}
		
		// Compute QY=QY-QX*(QX^T*QY)
		DenseMatrix Z=new DenseMatrix(QX.numRows(),QY.numColumns());
		DenseMatrix Temp1=new DenseMatrix(QX.numColumns(),QX.numRows());
		DenseMatrix Temp2=new DenseMatrix(QX.numColumns(),QY.numColumns());
		QX.transpose(Temp1); 	// Temp1=QX^T
		Temp1.mult(QY, Temp2);	// Temp2=Temp1*QY=QX^T*QY
		QX.mult(Temp2,Z);
		QY.add(-1D,Z);
		
		// Now compute the singular values
		try {
			SVD svd=SVD.factorize(QY);
			s=svd.getS();
		} catch (NotConvergedException e){}
		
		// Return angles in radians
		for (int i=0;i<s.length;++i) s[i]=Math.asin(s[i]);
		return s;
	}
	
	// Get sines of principal angles
	public static double[] getPrincipleAnglesSines(DenseMatrix X,DenseMatrix Y){
		double[] s=getPrincipleAngles(X,Y);
		for (int i=0;i<s.length;++i) s[i]=Math.sin(s[i]);
		return s;
	}
	
	// Get cosines of principal angles
	public static double[] getPrincipleAnglesCosines(DenseMatrix X,DenseMatrix Y){
		double[] s=getPrincipleAngles(X,Y);
		for (int i=0;i<s.length;++i) s[i]=Math.cos(s[i]);
		return s;
	}
	
	// Largest angle is in radians
	public static double getPrincipleAnglesLargest(DenseMatrix X,DenseMatrix Y){
		double[] s=getPrincipleAngles(X,Y);
		double max = 0;
		for (double d : s) if (d>max) max=d;
		return max;
	}
		

}




