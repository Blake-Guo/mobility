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

package sparse.eigensolvers.java;

import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.Matrices;
import org.junit.Test;
import sparse.eigenvolvers.java.Subspace;

public class SubspaceTest  extends TestCase {

	// Test accuracy for general angles
	@Test
	public void testAccuracyGeneral() {
		int n=100;
		// Input is tangent of angles
		double[] inputTangent={100,50,10,1,.5,1e-1,1e-4,1e-6,1e-8}; 
		
		DenseMatrix X=new DenseMatrix(n, inputTangent.length);
		DenseMatrix Y=new DenseMatrix(n, inputTangent.length);
		for (int i=0; i<inputTangent.length; ++i) X.set(i,i,1D);
		for (int i=0; i<inputTangent.length; ++i) X.set(i+inputTangent.length,i,inputTangent[i]);
		for (int i=0; i<inputTangent.length; ++i) Y.set(i,i,1D);
		
		// Multiply on right by random matrices (does not change column space)
		DenseMatrix XX=(DenseMatrix) Matrices.random(inputTangent.length,inputTangent.length);
		DenseMatrix YY=(DenseMatrix) Matrices.random(inputTangent.length,inputTangent.length);
		DenseMatrix Temp=new DenseMatrix(n, inputTangent.length);
		X.mult(XX, Temp); X=Temp;
		Temp=new DenseMatrix(n, inputTangent.length);
		Y.mult(YY, Temp); Y=Temp;
		
		double[] s=Subspace.getPrincipleAngles(X, Y);
		double err=0;
		for (int i=0; i<s.length; ++i) err=err+Math.pow(Math.atan(inputTangent[i])-s[i],2);
		err=Math.sqrt(err);
		assertTrue(err<1e-10);
	}
	
	
	// Test accuracy for small angles
	@Test
	public void testAccuracySmall() {
		int n=100;
		// Input is tangent of angles
		double[] inputTangent={1e-8,1e-16,1e-20,1e-30,0,0}; 
		
		DenseMatrix X=new DenseMatrix(n, inputTangent.length);
		DenseMatrix Y=new DenseMatrix(n, inputTangent.length);
		for (int i=0; i<inputTangent.length; ++i) X.set(i,i,1D);
		for (int i=0; i<inputTangent.length; ++i) X.set(i+inputTangent.length,i,inputTangent[i]);
		for (int i=0; i<inputTangent.length; ++i) Y.set(i,i,1D);
		
		double[] s=Subspace.getPrincipleAngles(X, Y);
		double err=0;
		for (int i=0; i<s.length; ++i) err=err+Math.pow(Math.atan(inputTangent[i])-s[i],2);
		err=Math.sqrt(err);
		assertTrue(err<1e-20);
	}
	
}
