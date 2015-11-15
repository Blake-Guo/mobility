package bike;

import java.io.*;
import java.util.*;

import no.uib.cipr.matrix.*;
import no.uib.cipr.matrix.sparse.ArpackSym;

import utility.*;

/*
 * Instruction for using this class
 * 
 * 	1: generate the laplacian matrix:	lap.genLaplacianMat_Dense("data/2014-07-trip data.csv");
 *  2: use the spectral theorem to reduce the dimension:		double[][] nvec = lap.spectralDimReduction(3);
 *	3: output the reduced matrix to a file:	lap.outputMatrix(nvec, "bike_vec_3d.csv", "\t");
 */

public class LaplacianGen {
	ODGraph odgraph;// the station map and od distribution
	public double[][] sim_mat;// the similarity matrix
	public double[] dia_mat;//diagonal matrix D
	public double[][] lap_mat;//laplacian matrix
	
	
	/**
	 * Generate the laplacian matrix.
	 * @param tripFile
	 * @param STFlag, to decide whether use temporal feature in the OD matrix or not.
	 * @throws IOException
	 */
	public void genLaplacianMat_Dense(String tripFile, boolean STFlag)throws IOException{
		
		//0: initialize the similairty matrix
		odgraph = new ODGraph(tripFile,STFlag);
		initSimMatrix_Dense();
		
		//1: construct the diagonal matrix and the Laplacian matrix
		dia_mat = new double[odgraph.stations.size()];
		lap_mat = new double[odgraph.stations.size()][];
		for(int i=0;i<lap_mat.length;i++)
			lap_mat[i] = new double[lap_mat.length];
		
		
		//2: compute the Diagonal matrix
		for(int i=0;i<dia_mat.length;i++)
		{
			dia_mat[i] = 0.0f;
			for(int j=0;j<dia_mat.length;j++)
				dia_mat[i] += sim_mat[i][j];
		}
		
		//3: compute the Laplacian matrix
		for(int i=0;i<lap_mat.length;i++){
			for(int j=i+1;j<lap_mat.length;j++)
				lap_mat[i][j] = lap_mat[j][i] = sim_mat[i][j] / Math.sqrt(dia_mat[i]) / Math.sqrt(dia_mat[j]);
		}
		
		System.out.println("Finish computing the laplacian matrix");
		
		return;
	}
	

	
	/**
	 * Generate the similarity matrix based on the OD_graph
	 */
	public void initSimMatrix_Dense() {

		// 1: construct the sim_mat
		sim_mat = new double[odgraph.stations.size()][];
		for (int i = 0; i < sim_mat.length; i++) {
			sim_mat[i] = new double[sim_mat.length];
		}

		// 2: compute each pair's similarity
		for(int i=0;i<sim_mat.length;i++)
			for(int j=i+1;j<sim_mat.length;j++)
				sim_mat[i][j] = sim_mat[j][i] = similarity(i,j);
		
		return;
	}

	
	/**
	 * Compute the similarity between two stations.
	 * @param sInd1
	 * @param sInd2
	 * @return
	 */
	public double similarity(int sInd1, int sInd2) {

		if (isNeighbor(sInd1, sInd2) == false)
			return 0.0;

		double cos_alpha = 0.6;
		double similarity = cos_alpha * cosineDistance(odgraph.od_mat[sInd1], odgraph.od_mat[sInd2]);
		similarity += (1-cos_alpha) * spatialSimilarity(sInd1,sInd2);

		return similarity;
	}
	
	
	/**
	 * compute the cosine distance between two vectors.
	 * @param a
	 * @param b
	 * @return
	 */
	public double cosineDistance(int[] a, int[] b){
		double ab = dot(a,b);
		double aa = dot(a,a);
		double bb = dot(b,b);
		
		return ab / (Math.sqrt(aa) * Math.sqrt(bb));
	}

	
	/**
	 * Compute the dot product of two arrays, they have must the same length;
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public double dot(int[] a, int[] b) {

		try {
			if (a.length != b.length)
				throw new Exception("two arrays have different lengths");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int sum = 0;
		for(int i=0;i<a.length;i++)
			sum += a[i] * b[i];

		return (double) sum;
	}

	
	/**
	 * Compute the spatial similarity value of two bike stations(given their indices, not id).
	 * 
	 * @param sInd1
	 * @param sInd2
	 * @return
	 */
	public double spatialSimilarity(int sInd1, int sInd2) {
		
		int id1 = odgraph.ind2id.get(sInd1);
		int id2 = odgraph.ind2id.get(sInd2);
		
		double dis = odgraph.stations.get(id1).distance_man(odgraph.stations.get(id2));
		
		double h = 600;//3-4 times of the average shortest distance
		
		return Math.exp(-dis/h);
	}
	
	

	/**
	 * Judge whether two stations sid1 and sid2 are neighbors(given their indices, not id).
	 * 
	 * @param sInd1
	 * @param sInd2
	 * @return
	 */
	public boolean isNeighbor(int sInd1, int sInd2) {
		double threshold = 2000;
		
		int id1 = odgraph.ind2id.get(sInd1);
		Station s1 = odgraph.stations.get(id1);
		
		int id2 = odgraph.ind2id.get(sInd2);
		Station s2 = odgraph.stations.get(id2);
		
		if(s1.distance_man(s2) < threshold)
			return true;
		
		return false;
	}
	
	
	/**
	 * Use the spectral theory to reduce the dimension to k and return the new data.
	 * Each row of matrix is the new vector 
	 * @param k
	 * @return
	 */
	public double[][] spectralDimReduction(int k){
		
		double[][] new_vecs = new double[sim_mat.length][];
		for(int i = 0;i<new_vecs.length;i++)
			new_vecs[i] = new double[k];
		
		Matrix lap = new DenseMatrix(lap_mat);
	
        ArpackSym solver = new ArpackSym(lap);
        
        Map<Double, DenseVectorSub> results = solver.solve(k, ArpackSym.Ritz.LA);
        
        int c=0;
        for (Map.Entry<Double, DenseVectorSub> e : results.entrySet()) {
        	for(int i=0;i<new_vecs.length;i++)
        		new_vecs[i][c] = e.getValue().get(i);
        	System.out.println(e.getKey());		
        	c++;
        }
        
        for(int i=0;i<new_vecs.length;i++)
        	Utilities.norm(new_vecs[i]);
        
        return new_vecs;
	}
	
	
	
	/**
	 * Compute the mobility similarity vs. distance.
	 * @param bin
	 * @return
	 */
	public double[] odSimvsDistance(int bin){
		int MAX_DIS = 10300;
		double[] cos_sim = new double[MAX_DIS / bin];
		int[] c = new int[MAX_DIS / bin];

		
		for(int i=0;i<odgraph.od_mat.length;i++)
		{
			for(int j=i+1;j<odgraph.od_mat.length;j++){
				double sim = cosineDistance(odgraph.od_mat[i],odgraph.od_mat[j]);
				Station ss = odgraph.stations.get( odgraph.ind2id.get(i));
				Station es = odgraph.stations.get( odgraph.ind2id.get(j));
				double dis = ss.distance_man(es);
				
				int ind = (int) dis/bin;
				if(ind >= c.length)
					ind = c.length - 1;
				
				cos_sim[ind] += sim;
				c[ind]++;
			}
		}
		
		
		for(int i=0;i<c.length;i++)
			cos_sim[i] /= c[i];
		
		
		return cos_sim;
	}
	
	/**
	 * For test, output the matrix to the specific file
	 * @param mat
	 * @param path
	 * @param sep, seperator '\t', ',' or ' '
	 * @throws IOException
	 */
	public void outputMatrix(double[][] mat, String path, String sep) throws IOException{
		BufferedWriter br = new BufferedWriter(new FileWriter(path));
		
		for(int i=0;i<mat.length;i++)
		{
			String line = Double.toString(mat[i][0]);
			
			for(int j=1;j<mat[i].length;j++)
				line += (sep + mat[i][j]);
			
			br.write(line);
			br.newLine();
		}
		
		br.close();
		
		System.out.println("Finish output the new feature vectors");
	}
	
	
	public static void main(String[] args)throws IOException{
		LaplacianGen lap = new LaplacianGen();
		lap.genLaplacianMat_Dense("data/2014-07-trip data.csv", true);
		
		//double[][] nvec = lap.spectralDimReduction(7);
		//lap.outputMatrix(nvec, "bike_vec.csv", "\t");
		
		
		
		double[] dist = lap.odSimvsDistance(150);
		for(int i=0;i<dist.length;i++)
			System.out.print(dist[i] + ",");
	}
}
