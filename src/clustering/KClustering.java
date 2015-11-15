package clustering;

import java.io.*;
import java.util.*;

import net.sf.javaml.clustering.Clusterer;
import net.sf.javaml.clustering.KMeans;
import net.sf.javaml.clustering.evaluation.ClusterEvaluation;
import net.sf.javaml.clustering.evaluation.SumOfSquaredErrors;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.Instance;
import net.sf.javaml.tools.data.FileHandler;

public class KClustering {

	public void outputClusters(String edgeVecsFile, String edgeInFile,
			String edgeClusOutFile, int k) throws Exception {
		Dataset data = FileHandler.loadDataset(new File("eigen_vecs.csv"));

		Clusterer km = new KMeans(30);

		Dataset[] clusters = km.cluster(data);

		
		/* Create a measure for the cluster quality */
		ClusterEvaluation sse = new SumOfSquaredErrors();

		/* Measure the quality of the clustering */
		double score = sse.score(clusters);
		
		
		System.out.println(k + ":" + score);
		
		System.out.println("Finish KNN Clustering");
		Map<Integer, Integer> cidMap = new HashMap<Integer, Integer>();

		for (int i = 0; i < clusters.length; i++) {
			for (int j = 0; j < clusters[i].size(); j++) {
				cidMap.put(clusters[i].get(j).getID(), i);
			}
		}

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
			int cid = cidMap.get(eid);
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

	public int bestClusterK(String edgeVecsInFile) throws Exception {

		Dataset data = FileHandler.loadDataset(new File(edgeVecsInFile));

		double minError = Double.MAX_VALUE;
		int bestK = 0;

		for (int i = 100; i <= 1000; i+=10) {
			Clusterer km = new KMeans(i);

			/* We cluster the data */
			Dataset[] clusters = km.cluster(data);

			/* Create a measure for the cluster quality */
			ClusterEvaluation sse = new SumOfSquaredErrors();

			/* Measure the quality of the clustering */
			double score = sse.score(clusters);
			
			
			System.out.println(i + ":" + score);
			
			if (score < minError) {
				minError = score;
				bestK = i;
			}

		}

		return bestK;
	}

	public static void main(String[] args) throws Exception {
		// Dataset data = FileHandler.loadDataset(new File("eigen_vecs.csv"));

		// Clusterer km = new KMeans(30);
		//
		// Dataset[] clusters = km.cluster(data);
		// for(int i=0;i<clusters.length;i++)
		// {
		// System.out.println(clusters[i].size() + ":");
		// for(int j=0;j<clusters[i].size() && j<5; j++)
		// {
		// System.out.print(clusters[i].get(j).getID()+":");
		// System.out.print(clusters[i].get(j)+"\t");
		// }
		// System.out.println();
		// }

		KClustering knn = new KClustering();
		
//		System.out.println(knn.bestClusterK("eigen_vecs.csv" ));
		knn.outputClusters("eigen_vecs.csv", "data/cedges.csv",
				"data/cedges_clustered.csv", 30);

	}

}
