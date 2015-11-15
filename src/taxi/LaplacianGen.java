package taxi;

import java.util.*;
import java.io.*;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVectorSub;
import no.uib.cipr.matrix.sparse.ArpackSym;
import no.uib.cipr.matrix.sparse.CompRowMatrix;
import road.*;
import sparse.eigenvolvers.java.Lobpcg;
import utility.*;

public class LaplacianGen {
	public Map<Integer, Double>[] sim_sparse_mat;
	public Map<Integer, Double>[] lap_sparse_mat;

	public ODGraph mgraph;

	/**
	 * Generate the similarity matrix based on the given mobility matrix file
	 * 
	 * @param nodeInFile
	 * @param edgeInFile
	 * @param mobSparseMatFile
	 * @throws IOException
	 */
	public void initSimMatrix_Sparse(String nodeInFile, String edgeInFile,
			String mobSparseMatFile) throws IOException {
		mgraph = new ODGraph();
		mgraph.readInRoadAndMGraph(nodeInFile, edgeInFile, mobSparseMatFile);

		double ave_cosim = 0;
		int icosim = 0;

		sim_sparse_mat = (HashMap<Integer, Double>[]) (new HashMap<?, ?>[mgraph.rnetwork.edges
				.size()]);

		for (int i = 0; i < mgraph.rnetwork.edges.size(); i++) {
			Set<Integer> neighs = mgraph.rnetwork.getEdgeKNeighbor(2, i);
			sim_sparse_mat[i] = new HashMap<Integer, Double>();

			for (int neigh : neighs) {
				double spatial_alpha = 0.05;// Maybe we will use a flexible
											// function to replace the constant.
				double cosim = dot_sparse(i, neigh) + spatial_alpha;

				if (neigh < i) {
					//to make sure the sparse matrix is symmetric, we choose the larger similarity value and assign it to both(means at least one of they two are the k neighbor of the other).
					double precosim = (sim_sparse_mat[neigh].containsKey(i) == false) ? 0.0
							: sim_sparse_mat[neigh].get(i);
					cosim = (precosim > cosim) ? precosim : cosim;
					sim_sparse_mat[neigh].put(i, cosim);
				}

				sim_sparse_mat[i].put(neigh, cosim);

				ave_cosim += cosim;
				icosim++;
			}
		}

		System.out.println("Finish generating the similarity matrix");
		System.out.println("The average similarity between each edge:"
				+ (ave_cosim / icosim));
	}

	/**
	 * Write the similarity matrix to a file.
	 * 
	 * @param outSimMatfile
	 * @throws IOException
	 */
	public void outputSimMatrix_Sparse(String outSimMatFile) throws IOException {

		BufferedWriter bw = new BufferedWriter(new FileWriter(outSimMatFile));

		int i = 0;
		System.out.println("edge number:" + mgraph.rnetwork.edges.size());
		while (i < mgraph.rnetwork.edges.size()) {// scan every edge
			StringBuilder strb = new StringBuilder();
			strb.append(i);

			for (Map.Entry<Integer, Double> entry : sim_sparse_mat[i]
					.entrySet()) {
				strb.append("\t" + entry.getKey() + ":" + entry.getValue());
				// strb.append('\t');
			}

			// strb.deleteCharAt(strb.length()-1);
			// System.out.println(strb.toString());
			bw.write(strb.toString());
			bw.newLine();
			bw.flush();

			i++;
		}

		bw.close();

		System.out.println("Finish outputing the sparse similarity matrix");

		return;
	}

	/**
	 * Calculate the cosine similarity between two edges
	 * 
	 * @param eid0
	 * @param eid1
	 * @return
	 */
	public double dot_sparse(int eid0, int eid1) {
		double sum = 0;
		Map<Integer, Integer> sparseVec0 = mgraph.m_sparse_mat[eid0];
		Map<Integer, Integer> sparseVec1 = mgraph.m_sparse_mat[eid1];

		for (Map.Entry<Integer, Integer> entry : sparseVec0.entrySet()) {
			int val1 = 0;
			if (sparseVec1.containsKey(entry.getKey()))
				val1 = sparseVec1.get(entry.getKey());
			sum += entry.getValue() * val1;
		}

		double len0 = 0;
		for (Map.Entry<Integer, Integer> entry : sparseVec0.entrySet()) {
			len0 += Math.pow(entry.getValue(), 2);
		}
		len0 = Math.sqrt(len0);

		double len1 = 0;
		for (Map.Entry<Integer, Integer> entry : sparseVec1.entrySet()) {
			len1 += Math.pow(entry.getValue(), 2);
		}
		len1 = Math.sqrt(len1);

		if (len0 == 0 || len1 == 0)
			return 0;
		return sum / (len0 * len1);
	}

	/**
	 * Generate the Laplacian matrix based on the given similarity matrix.
	 * 
	 * @param simSparseMatFile
	 * @throws IOException
	 */
	public void genLaplacianMat_Sparse(String simSparseMatFile)
			throws IOException {

		int nEdge = Utilities.countNumberOfLine(simSparseMatFile);
		sim_sparse_mat = (HashMap<Integer, Double>[]) (new HashMap<?, ?>[nEdge]);
		lap_sparse_mat = (HashMap<Integer, Double>[]) (new HashMap<?, ?>[nEdge]);
		double[] d = new double[nEdge];

		// 1: Read in the similarity matrix and generate the diagonal Matrix d.
		String line = "";
		BufferedReader br = new BufferedReader(new FileReader(simSparseMatFile));
		int eid = 0;
		while ((line = br.readLine()) != null) {

			String[] strs = line.split("\t");

			sim_sparse_mat[eid] = new HashMap<Integer, Double>();
			lap_sparse_mat[eid] = new HashMap<Integer, Double>();

			d[eid] = 0;
			for (int i = 1; i < strs.length; i++) {
				String[] keyval = strs[i].split(":");

				int tid = Integer.valueOf(keyval[0]);
				double tval = Double.valueOf(keyval[1]);
				d[eid] += tval;
				sim_sparse_mat[eid].put(tid, tval);
			}

			eid++;
		}
		br.close();

		// 2: Generate the laplacian matrix based on the similarity matrix and
		// diagonal matrix d.
		for (int i = 0; i < nEdge; i++) {

			// 1: D^(-1/2)*W*D^(-1/2)
			for (Map.Entry<Integer, Double> entry : sim_sparse_mat[i]
					.entrySet()) {

				int j = entry.getKey();
				double value = entry.getValue();

				// value = value / d[i];
				value = value / (Math.sqrt(d[i]) * Math.sqrt(d[j]));
				lap_sparse_mat[i].put(j, value);
			}

			lap_sparse_mat[i].put(i, 0.0);

			// //2: I - D^(-1/2)*W*D^(-1/2)
			// double dval = lap_sparse_mat[i].get(i);
			// dval = dval + 1.0;
			// lap_sparse_mat[i].put(i, dval);
		}

		System.out.println("Finish generating the Laplacian matrix");
	}

	/**
	 * Output the laplacian matrix to the file
	 * 
	 * @param outFile
	 * @throws IOException
	 */
	public void outputLaplacianMat_Sparse(String outFile) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));

		for (int i = 0; i < lap_sparse_mat.length; i++) {

			StringBuilder sb = new StringBuilder();
			for (Map.Entry<Integer, Double> entry : lap_sparse_mat[i]
					.entrySet()) {
				sb.append(entry.getKey() + ":" + entry.getValue() + "\t");
			}
			bw.write(sb.toString());
			bw.newLine();
			// bw.flush();
		}
		bw.close();

		System.out.println("Finish outputing the Laplacian matrix");
	}

	/**
	 * Read in the laplacian matrix
	 * 
	 * @param inFile
	 * @throws IOException
	 */
	public void readInLaplacianMat_Sparse(String inFile) throws IOException {

		BufferedReader br = new BufferedReader(new FileReader(inFile));
		int nEdge = Utilities.countNumberOfLine(inFile);
		lap_sparse_mat = (HashMap<Integer, Double>[]) (new HashMap<?, ?>[nEdge]);

		String line;
		int eid = 0;
		while ((line = br.readLine()) != null) {

			String[] strs = line.split("\t");

			lap_sparse_mat[eid] = new HashMap<Integer, Double>();

			for (int i = 0; i < strs.length; i++) {
				String[] keyval = strs[i].split(":");

				int tid = Integer.valueOf(keyval[0]);
				double tval = Double.valueOf(keyval[1]);
				lap_sparse_mat[eid].put(tid, tval);
			}
			eid++;
		}

		br.close();
	}

	public void eigenSolver_Lobpcg(String lapInFile, String eigenVecOutFile)
			throws IOException {
		readInLaplacianMat_Sparse(lapInFile);

		int nEdge = lap_sparse_mat.length;

		int[][] nz = new int[nEdge][];

		// fill nz
		for (int i = 0; i < nEdge; i++) {
			nz[i] = new int[lap_sparse_mat[i].size()];
			int j = 0;
			for (Map.Entry<Integer, Double> entry : lap_sparse_mat[i]
					.entrySet()) {
				nz[i][j++] = entry.getKey();
			}
			Arrays.sort(nz[i]);
		}

		CompRowMatrix crsMat = new CompRowMatrix(nEdge, nEdge, nz);

		for (int i = 0; i < nEdge; i++) {
			for (Map.Entry<Integer, Double> entry : lap_sparse_mat[i]
					.entrySet()) {
				crsMat.set(i, entry.getKey(), entry.getValue());
			}
		}

		Lobpcg lobpcg = new Lobpcg();
		lobpcg.setBlockSize(15);
		lobpcg.setMaxIterations(150);
		lobpcg.setResidualTolerance(5e-4);
		crsMat.scale(-1);

		// Run lobpcg with mostly default parameters
		lobpcg.runLobpcg(crsMat);

		// Get eigenvalues in double array
		System.out.println("Eigen values:");
		double[] eigenVals = lobpcg.getEigenvalues();
		for (int i = 0; i < eigenVals.length; i++) {
			System.out.print(eigenVals[i] + "\t");
		}

		// Output the eigen vectors
		DenseMatrix eigenvecs = lobpcg.getEigenvectors();
		BufferedWriter bw = new BufferedWriter(new FileWriter(eigenVecOutFile));
		for (int i = 0; i < eigenvecs.numRows(); i++) {
			StringBuilder sb = new StringBuilder();
			double[] rowVal = new double[eigenvecs.numColumns()];
			// Normalization
			double sum = 0;
			for (int j = 0; j < rowVal.length; j++) {
				rowVal[j] = eigenvecs.get(i, j);
				sum += Math.pow(rowVal[j], 2);
			}

			sum = Math.sqrt(sum);

			for (int j = 0; j < rowVal.length; j++) {
				rowVal[j] /= sum;
				sb.append(rowVal[j] + "\t");
			}
			bw.write(sb.toString());
			bw.newLine();
		}

		bw.close();
		// System.out.println("number of rows:" + eigenvecs.numRows());
		// System.out.println("number of cols:" + eigenvecs.numColumns());

		System.out.println("Finish the eigen solver");
	}

//	public void eigenSolver_MTJ(String lapInFile, String eigenVecOutFile)
//			throws IOException {
//		readInLaplacianMat_Sparse(lapInFile);
//
//		int nEdge = lap_sparse_mat.length;
//
//		int[][] nz = new int[nEdge][];
//
//		// fill nz
//		for (int i = 0; i < nEdge; i++) {
//			nz[i] = new int[lap_sparse_mat[i].size()];
//			int j = 0;
//			for (Map.Entry<Integer, Double> entry : lap_sparse_mat[i]
//					.entrySet()) {
//				nz[i][j++] = entry.getKey();
//			}
//			Arrays.sort(nz[i]);
//		}
//
//		CompRowMatrix crsMat = new CompRowMatrix(nEdge, nEdge, nz);
//
//		for (int i = 0; i < nEdge; i++) {
//			for (Map.Entry<Integer, Double> entry : lap_sparse_mat[i]
//					.entrySet()) {
//				crsMat.set(i, entry.getKey(), entry.getValue());
//			}
//		}
//
//		ArpackSym solver = new ArpackSym(crsMat);
//
//		Map<Double, DenseVectorSub> results = solver
//				.solve(2, ArpackSym.Ritz.LA);
//		
//		//Start to output the eigen value and vector
//		BufferedWriter bw = new BufferedWriter(new FileWriter(eigenVecOutFile));
//		for (Map.Entry<Double, DenseVectorSub> e : results.entrySet()) {
//
//			System.out.println(e.getKey());
//
//			DenseVectorSub evec = e.getValue();
//
//			// Output the eigen vectors
//
//			StringBuilder sb = new StringBuilder();
//			double[] rowVal = new double[evec.size()];
//			// Normalization
//			double sum = 0;
//			for (int i = 0; i < evec.size(); i++) {
//				rowVal[i] = evec.get(i);
//				sum += Math.pow(rowVal[i], 2);
//			}
//
//			sum = Math.sqrt(sum);
//
//			for (int j = 0; j < rowVal.length; j++) {
//				rowVal[j] /= sum;
//				sb.append(rowVal[j] + "\t");
//			}
//			bw.write(sb.toString());
//			bw.newLine();
//			bw.flush();
//		}
//
//		bw.close();
//		
//		System.out.println("Finish the eigen solver");
//	}

	public static void main(String[] args) throws IOException {
		LaplacianGen lg = new LaplacianGen();
		// lg.initSimMatrix_Sparse("data/cnodes.csv", "data/cedges.csv",
		// "sparse_mgraph.csv");
		// lg.outputSimMatrix_Sparse("sparse_sim_mat.csv");

		// lg.genLaplacianMat_Sparse("sparse_sim_mat.csv");
		// lg.outputLaplacianMat_Sparse("sparse_lap_mat.csv");

		//lg.eigenSolver_Lobpcg("sparse_lap_mat.csv", "eigen_vecs.csv");
	}
}
