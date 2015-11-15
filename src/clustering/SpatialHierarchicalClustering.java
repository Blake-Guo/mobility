package clustering;

import java.util.*;
import java.util.Map.Entry;

import road.*;

import java.io.*;

import no.uib.cipr.matrix.DenseVector;

//edge id and cluster id here start from 0.
public class SpatialHierarchicalClustering {
	Map<Integer, Set<Integer>> clusterSet;// for each cluster set, store the
											// edge ids belonging to it.
	Map<Integer, Integer> edge2ClusterIds;// given the edge id, return its
											// cluster id.
	RoadNetwork rnetwork;
	DenseVector[] eVecs;

	public SpatialHierarchicalClustering() {
	}

	public void initHierarchical(String nodeInFile, String edgeInFile,
			String edgeVecsFile) throws IOException {

		// 1: initialize the road network.
		rnetwork = new RoadNetwork(nodeInFile, edgeInFile);

		// 2: initialize the feature vector for each edge.
		eVecs = new DenseVector[rnetwork.edges.size()];
		initEdgeVecs(edgeVecsFile);

		// 3: Initialize the clustering coefficients
		clusterSet = new HashMap<Integer, Set<Integer>>();
		edge2ClusterIds = new HashMap<Integer, Integer>();

		for (int i = 0; i < rnetwork.edges.size(); i++) {
			edge2ClusterIds.put(i, i);

			Set<Integer> tset = new HashSet<Integer>();
			tset.add(i);
			clusterSet.put(i, tset);
		}

		System.out.println("Finish initializing the clustering coefficients");

	}

	/**
	 * Read in the edge vectors file and fill the mobility feature
	 * vector(DenseVector[] eVecs) for each edge. Before entering this function,
	 * we must make sure eVecs has been initialized through
	 * "eVecs = new DenseVector[rnetwork.edges.size()];".
	 * 
	 * @param edgeVecsFile
	 * @throws IOException
	 */
	void initEdgeVecs(String edgeVecsFile) throws IOException {

		BufferedReader br = new BufferedReader(new FileReader(edgeVecsFile));
		String line;
		int eid = 0;

		while ((line = br.readLine()) != null) {
			String[] strs = line.split("\t");

			eVecs[eid] = new DenseVector(strs.length);
			for (int i = 0; i < strs.length; i++)
				eVecs[eid].set(i, Double.valueOf(strs[i]));

			eid++;
		}

		br.close();

		System.out.println("Finish initializing and filling the edge vectors");
	}

	public void hierarchicalClustering(int nCluser) {
		double minDis;
		int id1 = -1, id2 = -1;
		double cThres = 5e-3;
		double cThres_Add = 1e-4;
		boolean fflag;

		while (clusterSet.size() > nCluser) {
			minDis = Double.MAX_VALUE;
			fflag = false;
			Set<Integer> vclusterIds = new HashSet<Integer>();

			for (Map.Entry<Integer, Set<Integer>> entry : clusterSet.entrySet()) {//go over all the clusters
				int cid1 = entry.getKey();
				vclusterIds.clear();

				for (Integer t_eid : entry.getValue()) {// go over each edge in
														// that cluster

					// go over each edge's neighborhood
					for (Integer neigh : rnetwork.edges.get(t_eid).adjEdgeIds) {

						int tcid2 = edge2ClusterIds.get(neigh);

						if (cid1 < tcid2 && vclusterIds.contains(tcid2) == false) {// we only cluster adjacent clusters
											// and avoid duplicate computation.
							vclusterIds.add(tcid2);
							double tdis = completeLinkage(cid1, tcid2);

							if (tdis < minDis) {// find a better cluster pair,
												// store their value.
								minDis = tdis;
								id1 = cid1;
								id2 = tcid2;
								// System.out.println("dis:" + minDis + "," +
								// id1 + ":" + id2);

								if (minDis < cThres) { // it is small enough,
														// merge them.
									fflag = true;
									break;
								}
							}

						}// if( cid1 < tcid2)
					}// for(Integer neigh : rnetwork.edges.get(t_eid).adjEdgeIds)

					if (fflag == true)
						break;

				}// go over each edge in that cluster
				if (fflag == true)
					break;
			}//go over each cluster

			// mergeTwoClusters, move the edges in id1 to id2.
			mergeTwoClusters(id1, id2);

			if (clusterSet.size() % 1000 == 0) {
				cThres += cThres_Add;
				if (cThres > 5e-2)
					cThres = 5e-2;
			}

			if (clusterSet.size() % 500 == 0) {
				System.out.println("cluster size:" + clusterSet.size());
			}

		}// while(clustersSet.size() > 30)

		System.out.println("Finish hierarchical clustering");
	}

	// public void hierarchicalClustering(int nCluser){
	// double minDis;
	// int id1=-1,id2=-1;
	// double cThres = 5e-3;
	// double cThres_Add = 1e-4;
	// boolean fflag;
	//
	// while(clusterSet.size() > nCluser){
	// minDis = Double.MAX_VALUE;
	// fflag = false;
	//
	// Iterator<Entry<Integer,Set<Integer>>> iter1 =
	// clusterSet.entrySet().iterator();
	//
	// while (iter1.hasNext()){
	// Entry<Integer,Set<Integer>> entry1 = iter1.next();
	// int cid1 = entry1.getKey();
	//
	// Iterator<Entry<Integer,Set<Integer>>> iter2 =
	// clusterSet.entrySet().iterator();
	// while(iter2.hasNext()){
	// Entry<Integer,Set<Integer>> entry2 = iter2.next();
	//
	// int tcid2 = entry2.getKey();
	//
	// if( cid1 < tcid2 && adjacent_cluster(cid1,tcid2)){//we only cluster
	// adjacent clusters and void duplicate computation.
	// double tdis = completeLinkage(cid1,tcid2);
	//
	// if(tdis < minDis){//find a better cluster pair, store their value.
	// minDis = tdis;
	// id1 = cid1;
	// id2 = tcid2;
	// //System.out.println("dis:" + minDis + "," + id1 + ":" + id2);
	//
	// if(minDis < cThres){ //it is small enough, merge them.
	// fflag = true;
	// break;
	// }
	// }
	//
	// }//if( cid1 < tcid2 && adjacent_cluster(cid1,tcid2))
	// }//while(iter2.hasNext())
	//
	// if(fflag == true)
	// break;
	//
	// }//while(iter1.hasNext())
	//
	// //mergeTwoClusters, move the edges in id1 to id2.
	// mergeTwoClusters(id1,id2);
	//
	//
	// if(clusterSet.size() % 1000 == 0)
	// {
	// cThres += cThres_Add;
	// if(cThres > 1e-2)
	// cThres = 1e-2;
	// }
	//
	// if(clusterSet.size() % 50 == 0){
	// System.out.println("cluster size:" + clusterSet.size());
	// }
	//
	// }//while(clustersSet.size() > 30)
	//
	//
	// System.out.println("Finish hierarchical clustering");
	// }

	/**
	 * Cluster two clusters, use cidi1 use the new merged cluster id and update
	 * clustersSet and edge2ClusterIds
	 * 
	 * @param cid1
	 * @param cid2
	 */
	public void mergeTwoClusters(int cid1, int cid2) {

		// 1: update the edge2ClusterIds
		Set<Integer> cluster2 = clusterSet.get(cid2);
		for (Integer eid : cluster2) {
			edge2ClusterIds.put(eid, cid1);
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

		for (Integer eid1 : cluster1)
			for (Integer eid2 : cluster2)
				if (adjacent_edge(eid1, eid2))
					return true;

		return false;
	}

	/**
	 * Given two edges, judge whether they are adjacent or not.
	 * 
	 * @param eid1
	 * @param eid2
	 * @return
	 */
	public boolean adjacent_edge(int eid1, int eid2) {

		Edge e1 = rnetwork.edges.get(eid1);
		for (Integer eid : e1.adjEdgeIds) {
			if (eid == eid2)
				return true;
		}

		return false;
	}

	/**
	 * Return the complete linkage between two clusters
	 * 
	 * @param cid1
	 * @param cid2
	 * @return
	 */
	public double completeLinkage(int cid1, int cid2) {

		double ldis = 0;

		Set<Integer> cluster1 = clusterSet.get(cid1);
		Set<Integer> cluster2 = clusterSet.get(cid2);

		for (Integer edgeId1 : cluster1) {
			for (Integer edgeId2 : cluster2) {
				double dis = eucDistance(eVecs[edgeId1], eVecs[edgeId2]);
				if (dis > ldis)
					ldis = dis;
			}
		}

		return ldis;
	}

	/**
	 * Return the Euclidean distance between two vectors v1 and v2.
	 * 
	 * @param v1
	 * @param v2
	 * @return
	 */
	public double eucDistance(DenseVector v1, DenseVector v2) {
		double sum = 0;
		for (int i = 0; i < v1.size(); i++) {
			sum += Math.pow(v1.get(i) - v2.get(i), 2);
		}

		sum = Math.sqrt(sum);
		return sum;
	}

	/**
	 * After we finish the clustering, we output the clustering result by adding
	 * the cluster id column to the original edge file.
	 * 
	 * @param edgeInFile
	 * @param edgeClusOutFile
	 * @throws IOException
	 */
	public void outputClusters(String edgeInFile, String edgeClusOutFile)
			throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(edgeInFile));
		BufferedWriter bw = new BufferedWriter(new FileWriter(edgeClusOutFile));

		String rline, wline;

		// title
		rline = br.readLine();
		wline = "cluster_id," + rline;
		bw.write(wline);
		bw.newLine();

		int eid = 0;
		while ((rline = br.readLine()) != null) {
			int cid = edge2ClusterIds.get(eid);
			wline = Integer.toString(cid) + "," + rline;
			bw.write(wline);
			bw.newLine();
			bw.flush();

			eid++;
		}

		bw.close();
		br.close();

		System.out.println("Finish writing the new edge clustering file");
	}

	public static void main(String[] args) throws IOException {
		SpatialHierarchicalClustering shcluster = new SpatialHierarchicalClustering();
		shcluster.initHierarchical("data/cnodes.csv", "data/cedges.csv",
				"eigen_vecs.csv");
		shcluster.hierarchicalClustering(30);
		shcluster.outputClusters("data/cnodes.csv",
				"data/cedges_sh_clustered.csv");

		// Map<Integer,Set<Integer>> map = new HashMap<Integer,Set<Integer>>();
		// map.put(5, new HashSet<Integer>());
		// map.put(8, new HashSet<Integer>());
		// map.remove(8);
		// Iterator<Integer> iter = map.keySet().iterator();
		// while(iter.hasNext())
		// System.out.println(iter.next());
	}
}
