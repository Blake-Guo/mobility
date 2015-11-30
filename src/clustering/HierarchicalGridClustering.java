package clustering;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.ArrayList;

import taxi_grid.Grid_Mobility;


public class HierarchicalGridClustering {
	Map<Integer, Set<Integer>> clusterSet;// for each cluster set, store the edge ids belonging to it.
	Map<Integer, Integer> grid2ClusterIds;// given the edge id, return its cluster id.
	ArrayList<ArrayList<Float>> feature_vecs;// the feature vectof of each neighborhood
	int grid_xnumber;
	int grid_ynumber;

	public HierarchicalGridClustering(String featureFile) throws IOException {
		clusterSet = new HashMap<Integer, Set<Integer>>();
		grid2ClusterIds = new HashMap<Integer, Integer>();
		feature_vecs = new ArrayList<ArrayList<Float>>();
		initHierarchical(featureFile);
	}

	/**
	 * Read in the feature vectors of each grid neighbor
	 * @param file
	 * @param sep
	 * @throws IOException
	 */
	private void readInGridFeatures(String file, String sep) throws IOException{
		BufferedReader breader = new BufferedReader(new FileReader(file));
		String line = "";

		while((line = breader.readLine()) != null){
			String[] strs = line.split(sep);
			ArrayList<Float> gridfeature = new ArrayList<Float>();
			for(String s:strs){
				gridfeature.add(Float.valueOf(s));
			}
			feature_vecs.add(gridfeature);
		}
		breader.close();
	}

	public void initHierarchical(String gridFeatureFile) throws IOException {

		// 1: initialize the grid feature vector.
		readInGridFeatures(gridFeatureFile, " ");
		Grid_Mobility grids = new Grid_Mobility();
		grid_xnumber = grids.grid_xnumber;
		grid_ynumber = grids.grid_ynumber;

		System.out.println("Grid X:"+grid_xnumber);
		System.out.println("Grid Y:"+grid_ynumber);

		// 2: Initialize the clustering
		for (int i = 0; i < feature_vecs.size(); i++) {
			grid2ClusterIds.put(i, i);

			Set<Integer> tset = new HashSet<Integer>();
			tset.add(i);
			clusterSet.put(i, tset);
		}
		System.out.println("Finish initializing the clustering coefficients");
	}

	/**
	 * Return the indices of the neibor grids of the given grid
	 * @param grid_ind
	 * @return
	 */
	public ArrayList<Integer> getAdjacentGrids(int grid_ind){
		ArrayList<Integer> neighbors = new ArrayList<Integer>();

		int x_ind = grid_ind % grid_xnumber; // --->x
		int y_ind = grid_ind / grid_xnumber; // | y

		int[] xoff = new int[]{0,0,1,-1};
		int[] yoff = new int[]{1,-1,0,0};

		for(int k=0;k<4;k++){
			int tmp_x = x_ind + xoff[k];
			int tmp_y = y_ind + yoff[k];

			if(tmp_x < 0 || tmp_x == grid_xnumber || tmp_y < 0 || tmp_y == grid_ynumber)
				continue;

			neighbors.add(tmp_y * grid_xnumber + tmp_x);
		}

		return neighbors;
	}

	/**
	 * The main clustering function.
	 * @param maximalClusterSize, the maximal number of grids in each cluster (to avoid a very larger cluster).
	 */
	public void hierarchicalClustering(int maximalClusterSize) {
		double minDis;
		int id1 = -1, id2 = -1;
		final double maximal_complete_dis = 1;
		//double cThres = 5e-3;
		//double cThres_Add = 1e-4;
		//boolean fflag;

		while (true) {
			minDis = Double.MAX_VALUE;//minimal complete distance between clusters
			//fflag = false;
			Set<Integer> vclusterIds = new HashSet<Integer>();
			id1 = -1;
			id2 = -1;

			for (Map.Entry<Integer, Set<Integer>> entry : clusterSet.entrySet()) {//go over all the clusters
				int cid1 = entry.getKey();
				vclusterIds.clear();

				for (Integer t_gid : entry.getValue()) {// go over each grid in
					// that cluster

					// go over each edge's neighborhood
					ArrayList<Integer> neighbors = getAdjacentGrids(t_gid);
					for (Integer neigh : neighbors) {

						int tcid2 = grid2ClusterIds.get(neigh);

						if(clusterSet.get(tcid2).size() + clusterSet.get(cid1).size() > maximalClusterSize)//these two clusters are too larger, we don't want to merge them
							continue;

						if (cid1 < tcid2 && vclusterIds.contains(tcid2) == false) {// we only cluster adjacent clusters
							// and avoid duplicate computation.
							vclusterIds.add(tcid2);
							double tdis = completeLinkage(cid1, tcid2);
							if (tdis < minDis) {// find a better cluster pair,
								// store their value.
								minDis = tdis;
								id1 = cid1;
								id2 = tcid2;
							}

						}// if( cid1 < tcid2)
					}// for(Integer neigh : rnetwork.edges.get(t_eid).adjEdgeIds)
				}// go over each edge in that cluster
			}//go over each cluster

			if(id1 == -1 || minDis > maximal_complete_dis)
				break;
			System.out.println("Minimal Complete Linkage:" + minDis );
			mergeTwoClusters(id1, id2);
		}
		System.out.println("Finish hierarchical clustering");
	}





	/**
	 * Cluster two clusters, use cidi1 use the new merged cluster id and update
	 * clustersSet and edge2ClusterIds
	 * 
	 * @param cid1
	 * @param cid2
	 */
	public void mergeTwoClusters(int cid1, int cid2) {
		// 1: update the grid2ClusterIds
		Set<Integer> cluster2 = clusterSet.get(cid2);
		for (Integer gid : cluster2) {
			grid2ClusterIds.put(gid, cid1);
		}

		// 2: update the clusterSet
		clusterSet.get(cid1).addAll(clusterSet.get(cid2));
		clusterSet.remove(cid2);
	}

	/**
	 * Given two clusters, judge whether they are adjacent or not.
	 * 
	 * @param cid1
	 * @param cid2
	 * @return
	 */
	public boolean adjacent_cluster(int cid1, int cid2) {
		Set<Integer> cluster1 = clusterSet.get(cid1);
		Set<Integer> cluster2 = clusterSet.get(cid2);

		for (Integer gid1 : cluster1)
			for (Integer gid2 : cluster2)
				if (adjacent_grid(gid1, gid2))
					return true;

		return false;
	}

	/**
	 * Given two grids, judge whether they are adjacent or not.
	 * 
	 * @param grid1
	 * @param grid2
	 * @return
	 */
	public boolean adjacent_grid(int grid1, int grid2) {

		int i1 = grid1 % grid_xnumber; // --->x
		int j1 = grid1 / grid_xnumber; // | y

		int i2 = grid2 % grid_xnumber; // --->x
		int j2 = grid2 / grid_xnumber; // | y

		if(Math.abs(i1-i2) + Math.abs(j1 - j2) == 1)
			return true;

		return false;
	}

	/**
	 * Return the complete linkage(farthest distance between two clusters) between two clusters.
	 * 
	 * @param cid1
	 * @param cid2
	 * @return
	 */
	public double completeLinkage(int cid1, int cid2) {
		double ldis = 0;
		Set<Integer> cluster1 = clusterSet.get(cid1);
		Set<Integer> cluster2 = clusterSet.get(cid2);

		for (Integer gId1 : cluster1) {
			for (Integer gId2 : cluster2) {
				double dis = cosDistance(feature_vecs.get(gId1), feature_vecs.get(gId2));
				if (dis > ldis)
					ldis = dis;
			}
		}
		return ldis;
	}


	/**
	 * Compute the cosine distance between two variables
	 * @param v1
	 * @param v2
	 * @return
	 */
	public double cosDistance(ArrayList<Float> v1, ArrayList<Float> v2){
		
		double cp = dotProduct(v1,v2);
		double squarev1 = dotProduct(v1,v1);
		double squarev2 = dotProduct(v2,v2);
		
		double dis = 1 - cp / (squarev1 * squarev2);
		
		return dis;
	}
	
	
	public double dotProduct(ArrayList<Float> v1, ArrayList<Float> v2){
		if(v1.size() != v2.size())
		{
			System.out.println("Error, two vectos have different size");
			return 0;
		}
		
		double sum = 0;
		for(int i=0;i<v1.size();i++)
			sum += v1.get(i) * v2.get(i);
		
		return sum;
	}


	public double eucDistance(ArrayList<Float> v1, ArrayList<Float> v2) {
		double sum = 0;
		for (int i = 0; i < v1.size(); i++) {
			sum += Math.pow(v1.get(i) - v2.get(i), 2);
		}
		sum = Math.sqrt(sum);
		return sum;
	}


	public void outputClusters(String clusterOutFile)
			throws IOException {

		BufferedWriter bw = new BufferedWriter(new FileWriter(clusterOutFile));

		String wline;
		int gridnumber = grid_xnumber * grid_ynumber;
		Map<Integer,Integer> clusterInd = new HashMap<Integer,Integer>();
		int ind = 0;
		for(int i=0;i<gridnumber;i++) {

			int cid = grid2ClusterIds.get(i);
			int tmp_ind = 0;
			if(clusterInd.containsKey(cid))
				tmp_ind = clusterInd.get(cid);
			else
			{
				clusterInd.put(cid, ind);
				tmp_ind = ind;
				ind++;
			}
			wline = Integer.toString(i) + " " + Integer.toString(tmp_ind);
			bw.write(wline);
			bw.newLine();
		}
		bw.close();
		System.out.println("Finish writing the new edge clustering file");
	}


	public static void main(String[] args) throws IOException {
		HierarchicalGridClustering gc = new HierarchicalGridClustering("../data/orig_features.txt");
		gc.hierarchicalClustering(6);
		gc.outputClusters("../data/cluster.txt");
	}

}
