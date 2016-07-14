
package Pathtime;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

import gnu.trove.map.hash.TShortDoubleHashMap;
import gnu.trove.map.hash.TShortShortHashMap;
import gnu.trove.iterator.TShortDoubleIterator;
import gnu.trove.set.hash.TShortHashSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.map.hash.TShortDoubleHashMap;
import gnu.trove.iterator.TShortIterator;
import gnu.trove.iterator.TShortDoubleIterator;
import gnu.trove.list.linked.TIntLinkedList;

import org.javatuples.Pair;
import org.javatuples.Quintet;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class Main {
	// Class to compare path based on their cost
	class PathCompare implements Comparator<Integer> {
		double[] pathCosts;
		public PathCompare(double[] pathCosts) {
			this.pathCosts = pathCosts;
		}

		public int compare(Integer p1, Integer p2) {
			double s1 = pathCosts[p1];
			double s2 = pathCosts[p2];
			// Descending order compare
			// As detailed here (https://docs.oracle.com/javase/7/docs/api/java/util/Comparator.html)
			// for ascending order, negative if first argument is less, 0 if 
			// equal and 1 if greater so the opposite for descending
			if (s1 < s2)
				return 1;
			else if (s1 > s2)
				return -1;
			else
				return 0;
		}
	}
	
	class GeneObjDeltaCompare implements Comparator<Integer> {
		double [] geneCosts;
		public GeneObjDeltaCompare(double[] geneCosts) {
			this.geneCosts = geneCosts;
		}
		
		public int compare(Integer g1, Integer g2) {
			double s1 = geneCosts[g1];
			double s2 = geneCosts[g2];
			
			if (s1 < s2)
				return 1;
			else if (s1 > s2)
				return -1;
			else
				return 0;
		}
	}

	public double l1, l2;
	public int L, N;
	public short[][] P; // paths
	public double[] w_P; // path scores
	public int[][] P_g;
	public byte[] b_p, b_g; // bv for paths, genes
	public byte[] best_b_p, best_b_g; // bv for paths, genes
	public short[] index_to_targets;
	public short[] targets_to_index;
	public short[][] T;
	public double[][] f_g_t;
	public double[] f;
	public double[] delta_g;

	public TShortDoubleHashMap[] edgeScores;
	public boolean[] isTf;
	public boolean[] isMirna;
	public boolean[] isRnaHit;
	public boolean [] isSource;
	public short [][] nl;
	public short[] sources;

	public int numGenes;

	public HashMap<String, Short> g_to_i; // gene name to index
	public String [] i_to_g; // gene index to name
	public String [] old_i_to_g; // gene index to name

	public int numPaths;
	public int numPhases;
	public TShortDoubleHashMap[][] tfToGene;
	public TShortHashSet[] phaseGenes;
	public int pathBufferSize;
	public int maxNumPaths;
	public int maxPathLength;
	public int maxDegree;
	public String mirnaPrefix;
	public String pathFile;
	public int numPhaseGenes;
	public double[][] foldChange;
	public TShortHashSet[] phaseTfs;

	public static final boolean DEBUG = true;
	
	public static HashMap<String, String> readConfigFile(String configFilename) {
		HashMap<String, String> args = new HashMap<String, String>();
		BufferedReader br = null;
		try {
			br = new BufferedReader(new FileReader(configFilename));
			String line;
			int i = 0;
			for (; (line = br.readLine()) != null;) {
				i += 1;
				if (line.trim() == "") {
					continue;
				}
				String [] lineSplit = line.split("#");
				String temp = lineSplit[0].trim();
				if (lineSplit.length >= 2 && temp.length() == 0) {
					continue;
				}
				if (lineSplit[0].trim().length() == 0) {
					continue;
				}
				lineSplit = lineSplit[0].trim().split("\t");
				if (lineSplit.length == 0) {
					continue;
				}
				if (lineSplit.length == 1) {
					System.err.println("ERROR: parameter " + lineSplit[0] + " in config file " + configFilename + " does not have a value! (line " + i + ")");
					System.exit(1);
				}
				
				System.err.println(lineSplit[0] + "\t" + lineSplit[1]);
				args.put(lineSplit[0].trim().toLowerCase(), lineSplit[1].trim());
			}
			br.close();
		}
		catch(IOException ie) {
			ie.printStackTrace();
		}

		return args;
	}

	public static void main(String[] args) {

		if (args.length < 1) {
			System.err.println("Usage: Main.java <config file>");
			System.exit(1);
		}
		
		HashMap<String, String> h = readConfigFile(args[0]);
		Main main = new Main(h);
	}
	
	public static String argValue(
			HashMap<String, String> args,
			String argName,
			String defaultValue,
			boolean mandatory)
	{
		String temp = argName;
		argName = argName.toLowerCase();
		if (args.containsKey(argName)) {
			return args.get(argName);
		}
		else if (mandatory) {
			throw new RuntimeException(temp + " is mandatory but not provided");
		}
		return defaultValue;
	}
	
	public void computeImportantGenes(HashMap<String, String> args) {
		double minRankingFoldChange = Double.parseDouble(
				argValue(args, "minRankingFoldChange", "2", false));
		String sourceFilepath = argValue(args, "sourceFilepath", "", true);
	}
	
	public Main(HashMap<String, String> args) {
		
		pathFile = argValue(args, "pathFile", "", true);
		numPhaseGenes = Integer.parseInt(argValue(args, "numTPtargets", "", true));
		maxNumPaths = Integer.parseInt(argValue(args, "maxNumPaths", "100000", false));
		maxPathLength = Integer.parseInt(argValue(args, "maxPathLength", "10", false));

		pathBufferSize = maxNumPaths;
		mirnaPrefix = argValue(args, "mirnaPrefix", "hsa", false);
		l1 = Double.parseDouble(argValue(args, "l1", "1", false));
		l2 = Double.parseDouble(argValue(args, "l2", "0.3", false));
		L = Integer.parseInt(argValue(args, "L", "200", false));
		N = Integer.parseInt(argValue(args, "N", "1000", true));

		System.err.println("pathFile: " + pathFile);

		P = new short[50000000][];
		w_P = new double[50000000];
		for (int i=0; i<w_P.length; ++i) {
			w_P[i] = -Double.MAX_VALUE;
		}
		numPaths = 0;

		System.err.println("Reading data");
		readData(args);
		System.err.println("Data read");

		searchPaths();

		// Sort searched paths
		Integer[] pathIndexList = new Integer[numPaths];
		for (int i=0; i<numPaths; ++i) {
			pathIndexList[i] = i;
		}
		Arrays.sort(pathIndexList, new PathCompare(w_P));
		double totalScore = 0;
		for (int i=0; i<numPaths; ++i) {
			// Swap paths
        		short[] temp = P[i];
        		double temp2 = w_P[i];
        		P[i] = P[pathIndexList[i]];
        		w_P[i] = w_P[pathIndexList[i]];
        		P[pathIndexList[i]] = temp;
        		w_P[pathIndexList[i]] = temp2;
		}

		// Add targets to paths
		addTargetsToPaths();

		for (int i=0; i<numPaths; ++i) {
			w_P[i] = Math.exp(w_P[i]);
		}

		// reasses numGenes and compute total path weight;
		g_to_i.clear();
		assert g_to_i.size() == 0 : g_to_i.size();
		numGenes = 0;
		for (int i=0; i<numPaths; ++i) {
			for (int k=0; k<P[i].length; ++k) {
				String gene = i_to_g[P[i][k]];
				if (!g_to_i.containsKey(gene)) {
					short geneIndex = (short) g_to_i.size();
					g_to_i.put(gene, geneIndex);
				}
				short geneindex = g_to_i.get(gene);
				P[i][k] = geneindex;
			}
		}

		numGenes = g_to_i.size();
		old_i_to_g = i_to_g;
		i_to_g = new String[numGenes];
		for (Map.Entry<String, Short> entry : g_to_i.entrySet()) {
			assert entry.getValue() < i_to_g.length : 
				entry.getValue() + " < " + i_to_g.length;
			i_to_g[entry.getValue()] = entry.getKey();
		}
		
		System.err.println("Final num genes: " + numGenes + "\n\n");

		TShortHashSet targetSet = new TShortHashSet();
		int[] numPathsPerGene = new int[numGenes];
		for (int i=0; i<numPaths; ++i) {
			totalScore += w_P[i];
			assert P[i] != null;
			assert targetSet != null;
			targetSet.add(P[i][P[i].length-1]);
			for (int k=0; k<P[i].length; ++k) {
				++numPathsPerGene[P[i][k]];
			}
		}
		System.err.println("totalScore: " + totalScore);
		
		l1 *= totalScore/numGenes;
		l2 *= totalScore/targetSet.size();
		delta_g = new double[numGenes];
		b_p = new byte[numPaths];
		best_b_p = new byte[numPaths];
		for (int i=0; i<numPaths; ++i) {
			b_p[i] = 1;
		}
		b_g = new byte[numGenes];
		best_b_g = new byte[numGenes];
		P_g = new int[numGenes][];
		int[] genePathIndex = new int[numGenes];
		targets_to_index = new short[numGenes];
		for (int i=0; i<numGenes; ++i) {
			b_g[i] = 1;
			P_g[i] = new int[numPathsPerGene[i]];
			targets_to_index[i] = -1;
			genePathIndex[i] = numPathsPerGene[i];
		}
		for (int i=0; i<numPaths; ++i) {
			for (int k=0; k<P[i].length; ++k) {
				int g = P[i][k];
				P_g[g][--genePathIndex[g]] = i;
			}
		}
		index_to_targets = new short[targetSet.size()];
		T = new short[index_to_targets.length][];

		short k=0;
		for (TShortIterator it=targetSet.iterator(); it.hasNext(); ++k) {
			short target = it.next();
			index_to_targets[k] = target;
			targets_to_index[target] = k;
		}

		ArrayList<HashSet<Short>> targets_to_genes = 
				new ArrayList<HashSet<Short>>();
		for (int i=0; i<index_to_targets.length; ++i) {
			targets_to_genes.add(new HashSet<Short>());
		}
		for (int i=0; i<numPaths; ++i) {
			short target = P[i][P[i].length-1];
			int targetIndex = targets_to_index[target];
			for (int j=0; j<P[i].length-1; ++j) {
				targets_to_genes.get(targetIndex).add(P[i][j]);
			}
			targets_to_genes.get(targetIndex).add(target);
		}

		f_g_t = new double[index_to_targets.length][];
		f = new double[index_to_targets.length];
		for (int i=0; i<index_to_targets.length; ++i) {
			T[i] = new short[targets_to_genes.get(i).size()];
			f_g_t[i] = new double[numGenes];
			int j = 0;
			for (short g : targets_to_genes.get(i)) {
				T[i][j++] = g; 
			}
		}

		boolean filter = argValue(args, "filterPaths?", "true", false).toLowerCase() == "true";
		if (filter) {
			filterPaths();
		}

		printPaths();
	}
	
	public void addTargetsToPaths() {
		int [] numPathsForSource = new int[numGenes];
		int [] numPathsForTarget = new int[numGenes];
		int [] numPathsForTf = new int[numGenes];
		int curP = P.length-1;

		int temp = 0;
		for (int i=0; i<numPaths; ++i) {
			short [] p = P[i];
			int l = p.length;
			Pair<Integer, Integer> phaseBounds = getStartEndPhase(p, l);
			int startPhase = phaseBounds.getValue0();
			int endPhase = phaseBounds.getValue1();

			short tf = p[p.length-1];
			short source = p[0];
			
			assert isSource[source];
			assert isTf[tf];

			//System.err.println("i: " + i + ", source: " + source 
			//		+ ", tf: " + tf + ", startPhase: " + startPhase 
			//		+ ", endPhase: " + endPhase);
			//temp = 0;
			for (int j=startPhase; j<=endPhase; ++j) {
				assert tfToGene[tf][j] != null;
				if (phaseTfs[j] != null && !phaseTfs[j].contains(tf)) {
					continue;
				}
				for (TShortDoubleIterator it=tfToGene[tf][j].iterator();
						it.hasNext();)
				{
					it.advance();
					short gene = it.key();
					
					/*
					if (isMirna[tf] 
							&& 1 == (int) (
									Math.signum(foldChange[tf][startPhase])
									*Math.signum(foldChange[gene][j])))
					{
						continue;
					}
					*/

					double score = it.value();
					if (((numPathsForSource[source] >= maxNumPaths) || (numPathsForTarget[gene] >= maxNumPaths)) && maxNumPaths != 0)
					{
						continue;
					}

					boolean use = true;
					for (int k=0; k<l; ++k) {
						if (p[k] == gene) {
							use = false;
							break;
						}
					}

					if (!use) {
						break;
					}


					//++temp;
					++numPathsForTarget[gene];
					++numPathsForSource[source];
					++numPathsForTf[tf];
					P[curP] = new short[l+1];
					w_P[curP] = w_P[i]+score;
					for (int k=0; k<l; ++k) {
						P[curP][k] = p[k];

					}
					P[curP][l] = gene;
					--curP;
				}
			}
			//System.err.println("source: " + source 
			//		+ ", numPathsForSource: " + temp);
		}
		numPaths = P.length-curP-1;
		System.err.println("Final numPaths: " + numPaths);
		
		for (int i=0; i<numPaths; ++i) {
			P[i] = P[P.length-i-1];
			w_P[i] = w_P[P.length-i-1];
		}
	}

	// Get start and end phase indices of the phases to be explored for
	// the target
	public Pair<Integer, Integer> getStartEndPhase(short [] p, int l) {
		int pathPhase = -1;
		int firstCurPhaseGene = 0;

		// pathPhase is set to the highest phase gene in the path
		// firstCurPhaseGene is the index of the earliest gene in the path
		// belonging to a phase <= pathPhase
		for (int k=l-1; k>=0; --k) {
			short gene = p[k];
			for (int j=numPhases-1; j>=0; --j) {
				if (phaseGenes[j].contains(gene)) {
					if (pathPhase < j) {
						pathPhase = j;
					}
					if (!(pathPhase > j)) {
						firstCurPhaseGene = k;
					}
				}
			}
		}

		// If previous phase gene present before EARLIEST gene of current
		// phase
		boolean hasPrevPhaseGene = false;
		if (pathPhase > 0) {
			for (int k=firstCurPhaseGene-1; k>=0; --k) {
				short gene = p[k];
				if (phaseGenes[pathPhase-1].contains(gene)) {
					hasPrevPhaseGene = true;
					break;
				}
			}
		}

		assert pathPhase < numPhases;
		assert pathPhase >= -1;
		int startPhase = pathPhase;
		int endPhase = pathPhase;

		if (pathPhase <= 0) {
			startPhase = 0;
			endPhase = 0;
		}
		else if (pathPhase > 0 && !hasPrevPhaseGene) {
			startPhase = pathPhase+1;
			endPhase = pathPhase+1;
		}
		else {
			endPhase = pathPhase+1;
		}

		endPhase = Math.min(endPhase, numPhases-1);

		return new Pair<Integer, Integer>(startPhase, endPhase);
	}

	public void searchPaths() {
		int actualPathCapacity = pathBufferSize*2 + maxDegree;

		int [][] pathLength = {new int[actualPathCapacity],
				new int[actualPathCapacity]};
		double [][] pathCosts = {new double[actualPathCapacity],
				new double[actualPathCapacity]};
		short [][][] pathBuffer = {
				new short[actualPathCapacity][maxPathLength],
				new short[actualPathCapacity][maxPathLength]};
		Integer [] pathIndexList = new Integer[actualPathCapacity];

		int[] numPathsForSource = new int[numGenes];
		int[] numPathsForTarget = new int[numGenes];                
		int numPathsInBuffer =
				setupPaths(
						pathIndexList,
						pathLength,
						pathCosts,
						pathBuffer);

		int numPathsDone = 0;

		System.err.println("Compute paths");
		int curP;
		for (; ;) {
			Quintet<Integer,
			Integer,
			Integer,
			Integer,
			Integer> ret1 = fillPathBufferAndWritePaths(
					numPathsInBuffer,
					numPathsDone,
					pathLength[0],
					pathBuffer[0],
					pathCosts[0],
					-1,
					numPathsForSource,
					numPathsForTarget);
			curP = ret1.getValue0();
			numPathsInBuffer = ret1.getValue2();
			numPathsDone = ret1.getValue3();
			int numNewBufferPaths = ret1.getValue4();
			System.err.println(numPathsDone + " paths done");
			System.err.println(numNewBufferPaths + " numNewBufferPaths");
			if (numNewBufferPaths == 0) {
				break;
			}

			int ret2 = fixPathBuffer(
					pathLength,
					pathCosts,
					pathBuffer,
					pathIndexList,
					numPathsInBuffer,
					curP,
					actualPathCapacity
					);
			numPathsInBuffer = ret2;
			if (numPathsInBuffer == 0) {
				break;
			}
		}
		System.err.println(numPathsDone + " total init search paths done");
		System.err.println("Search Done!");
	}

	int setupPaths(
			Integer [] pathIndexList,
			int [][] pathLength,
			double [][] pathCosts,
			short [][][] pathBuffer
			) {
		assert pathLength.length == 2 : pathLength.length;
		assert pathCosts.length == 2 : pathCosts.length;
		assert pathBuffer.length == 2 : pathBuffer.length;

		assert pathLength[0].length == pathLength[1].length;
		assert pathCosts[0].length == pathCosts[1].length;
		assert pathBuffer[0].length == pathBuffer[1].length;
		assert pathLength[0].length == pathIndexList.length;

		// Setup paths
		for (int i=0; i<pathCosts[0].length; ++i) {
			pathIndexList[i] = i;
			pathLength[0][i] = 0;
			pathLength[1][i] = 0;
			pathCosts[0][i] = -Double.MAX_VALUE;
			pathCosts[1][i] = -Double.MAX_VALUE;
		}

		int numSources = 0;
		for (short i=0; i<numGenes; ++i) {
			if (isTf[i] && nl[i].length > 0) {
				pathLength[0][numSources] = 1;
				pathCosts[0][numSources] = 0;
				pathBuffer[0][numSources++][0] = i;
			}
		}
		return numSources;
	}

	/*
    Sorts the pathbuffer so that the best paths are first and readies
    it for a new search
	 */
	int fixPathBuffer(
			int [][] pathLength,
			double [][] pathCosts,
			short [][][] pathBuffer,
			Integer [] pathIndexList,
			int numPathsInBuffer,
			final int curP,
			final int actualPathCapacity
			) {
		assert curP > 0;
		assert pathLength.length == 2 : pathLength.length;
		assert pathCosts.length == 2 : pathCosts.length;
		assert pathBuffer.length == 2 : pathBuffer.length;

		assert pathLength[0].length == pathLength[1].length;
		assert pathCosts[0].length == pathCosts[1].length;
		assert pathBuffer[0].length == pathBuffer[1].length;
		assert pathLength[0].length == pathIndexList.length;

		if (DEBUG) {
			for (int i=0; i<numPathsInBuffer; ++i) {
				assert pathLength[0][i] > 0
				: i + "\t" + numPathsInBuffer + "\t"+ curP;
				assert pathCosts[0][i] > -Double.MAX_VALUE;
			}
			for (int i=0; i<actualPathCapacity; ++i) {
				if (pathLength[0][i] == 0) {
					assert pathCosts[0][i] == -Double.MAX_VALUE : i;
				}
			}
		}

		// Everything upto (but not including curP as been processed into
		// new paths. So set the path costs of those paths to the minimum
		// so on sorting they come in last. And set their path lengths to 0
		for (int k=0; k<curP; ++k) {
			pathCosts[0][k] = -Double.MAX_VALUE;
			pathLength[0][k] = 0;
		}

		// Sort pathIndexList according to pathCosts[0]
		for (int i=0; i<actualPathCapacity; ++i) {
			pathIndexList[i] = i;
		}
		Arrays.sort(pathIndexList, new PathCompare(pathCosts[0]));

		// an element of pathIndexList is the index of the path in the old
		// path array. Thus if the value of an element is 46, it might be
		// in the 1000 position in pathIndexList but it refers to the the
		// 46th path in pathBuffer[0], pathCosts[0], and pathLength[0]

		if (DEBUG) {
			int i=1;
			for (i=1; i<actualPathCapacity; ++i) {
				assert (pathCosts[0][pathIndexList[i]]
						<= pathCosts[0][pathIndexList[i-1]]);
			}
			for (i=0; i<numPathsInBuffer-curP; ++i) {
				assert pathLength[0][pathIndexList[i]] > 0
				: i + "\t" + pathIndexList[i]
						+ "\t" + numPathsInBuffer + "\t"+ curP
						+ "\t" + pathCosts[0][pathIndexList[i]];
				assert pathCosts[0][pathIndexList[i]] > -Double.MAX_VALUE;
			}
		}

		// Reorder the pathBuffer[0], pathCosts[0], and pathLength[0] according
		// to pathIndexList
		for (int k=0; k<actualPathCapacity; ++k) {
			pathBuffer[1][k] = pathBuffer[0][pathIndexList[k]];
			pathCosts[1][k] = pathCosts[0][pathIndexList[k]];
			pathLength[1][k] = pathLength[0][pathIndexList[k]];
		}

		short[][] pathBufferTemp = pathBuffer[0];
		pathBuffer[0] = pathBuffer[1];
		pathBuffer[1] = pathBufferTemp;

		double[] pathCostTemp = pathCosts[0];
		pathCosts[0] = pathCosts[1];
		pathCosts[1] = pathCostTemp;

		int[] pathLengthTemp = pathLength[0];
		pathLength[0] = pathLength[1];
		pathLength[1] = pathLengthTemp;

		int newNumPathsInBuffer = Math.min(numPathsInBuffer-curP, pathBufferSize);

		return newNumPathsInBuffer;
	}
	Quintet<Integer, Integer, Integer, Integer, Integer>
	fillPathBufferAndWritePaths(
			int numPathsInBuffer,
			int numPathsDone,
			int[] pathLength,
			short[][] pathBuffer,
			double[] pathCosts,
			final int curGeneIndex
			) {
		return fillPathBufferAndWritePaths(
				numPathsInBuffer, 
				numPathsDone, 
				pathLength, 
				pathBuffer, 
				pathCosts, 
				curGeneIndex,
				null,
				null);
	}

	/*
    Searches for new paths until the search buffer fills up
	 */
	Quintet<Integer, Integer, Integer, Integer, Integer>
	fillPathBufferAndWritePaths(
			int numPathsInBuffer,
			int numPathsDone,
			int[] pathLength,
			short[][] pathBuffer,
			double[] pathCosts,
			final int curGeneIndex,
			int[] numPathsForSource,
			int[] numPathsForTarget
			) {
		if (numPathsForSource == null) {
			throw new RuntimeException("numPathsForSource cannot be null");
		}
		if (numPathsForTarget == null) {
			throw new RuntimeException("numPathsForTarget cannot be null");
		}

		// This loop tries expanding the path indexed by curP according
		// to the edge indexed by curE in the nl array of the last node
		// of curP
		//long curTime = System.currentTimeMillis();
		int newBufferPaths = 0;

		int curP = 0;
		int curE = 0;
		if (DEBUG) {
			for (int i=0; i<pathLength.length; ++i) {
				if (pathLength[i] == 0) {
					assert pathCosts[i] == -Double.MAX_VALUE : i;
				}
			}
		}
		int numTimes = 0;
		assert numPathsInBuffer < pathBufferSize*2;
		assert numPathsInBuffer > 0;
		for (curP = 0;
				curP < numPathsInBuffer
				&& numPathsInBuffer<pathBufferSize*2;
				++curP)
		{
			if (pathLength[curP] == 0) {
				if (DEBUG) {
					for (int i=curP+1; i<pathBufferSize*2; ++i) {
						assert pathLength[i] == 0;
					}
				}
				System.err.println("pathLength break");
				assert pathCosts[curP] == -Double.MAX_VALUE;
				++curP;
				break;
			}
			assert pathLength[curP] > 0;

			// We should never put full length paths into the buffer
			assert pathLength[curP] < maxPathLength;

			final short endGene = pathBuffer[curP][pathLength[curP]-1];
			assert nl[endGene].length > 0 : curP;
			assert pathCosts[curP] > -Double.MAX_VALUE;

			// At regular time intervals print numPathsDone
			//long curTime2 = System.currentTimeMillis();
			//if ((curTime2-curTime) > 10*1000) {
			//   System.err.println(numPathsDone + " paths done");
			//   curTime = curTime2;
			//}

			NewPathGenes:
				for (curE=0; curE < nl[endGene].length; ++curE) {
					if (curE > 0) {
						assert nl[endGene][curE] != nl[endGene][curE-1] :
							nl[endGene][curE]; 

					}
					final short newGene = nl[endGene][curE];
					for (int k=0; k<pathLength[curP]; ++k) {
						if (pathBuffer[curP][k] == newGene) {
							continue NewPathGenes;
						}
					}

					// If we reach here, that means adding newGene doesn't result
					// in a cycle.

					// Copy new path to the position numPathsInBuffer and set
					// new path length and path cost
					// But don't increment numPathsInBuffer yet as we may
					// not wanna add this path to the buffer if it's of maxlength
					for (int k=0; k<pathLength[curP]; ++k) {
						pathBuffer[numPathsInBuffer][k] = pathBuffer[curP][k];
					}
					pathBuffer[numPathsInBuffer][pathLength[curP]] = newGene;
					assert pathLength[curP] > 0;
					pathLength[numPathsInBuffer] = pathLength[curP]+1;
					pathCosts[numPathsInBuffer] =
							pathCosts[curP]
									+ edgeScores[endGene].get(newGene);
					assert pathCosts[numPathsInBuffer] > -Double.MAX_VALUE;

					int newPathIndex = numPathsInBuffer;

					// If new gene is a target, then print it.
					if (isSource[newGene]) {
						short tf = newGene;
						short source = pathBuffer[newPathIndex][0];

						if ((numPathsForSource[source] < maxNumPaths
								&& numPathsForTarget[tf] < maxNumPaths)
								&& maxNumPaths != 0)
						{
							short [] p = pathBuffer[newPathIndex];
							P[numPaths] = new short[pathLength[newPathIndex]];
							int l = pathLength[newPathIndex];
							for (int k=0; k<l; ++k) {
								P[numPaths][k] = p[l-k-1];
							}
							w_P[numPaths] = pathCosts[newPathIndex];
							++numPaths;
							++numPathsDone;
							++numPathsForTarget[tf];
							++numPathsForSource[source];
						}
						else {
							++numTimes;
						}
					}

					// If new path can be expanded, then the new path is 
					// registered in the buffer
					if (pathLength[newPathIndex] < maxPathLength
							&& nl[newGene].length > 0) {
						++numPathsInBuffer;
						++newBufferPaths;
					}
					else {
						pathLength[newPathIndex] = 0;
						pathCosts[newPathIndex] = -Double.MAX_VALUE;
					}
				}
		}
		assert curP > 0;
		if (numTimes > 0) {
			System.err.println("Num Times: " + numTimes);
		}
		if (DEBUG) {
			for (int i=0; i<pathLength.length; ++i) {
				if (pathLength[i] == 0) {
					assert pathCosts[i] == -Double.MAX_VALUE : i;
				}
			}
		}
		//System.err.println("curP " + curP + ", curE: " + curE);
		return new Quintet<Integer, Integer, Integer, Integer, Integer>(
				curP,
				curE,
				numPathsInBuffer,
				numPathsDone,
				newBufferPaths
				);
	}

	public void printPaths() {
		BufferedWriter br = null;
		try {
			br = new BufferedWriter(new FileWriter(pathFile));
			for (int i=0; i<numPaths; ++i) {
				String s = Double.toString(Math.exp(w_P[i]))
						+ "\t" + P[i].length;
				for (int k=P[i].length-1; k>=0; --k) {
					s += "\t" + i_to_g[P[i][k]];
				}
				if (b_p[i] == 1) {
					br.write(s + "\n");
				}
			}
			br.close();
		}
		catch (IOException ie) {
			ie.printStackTrace();
		}

		br = null;
		try {
			br = new BufferedWriter(new FileWriter(pathFile + ".phasegenes.txt"));

			TShortIterator [] its = new TShortIterator[phaseGenes.length];
			int size=0;
			for (int i=0; i<phaseGenes.length; i++) {
				size = Math.max(size, phaseGenes[i].size());
				its[i] = phaseGenes[i].iterator();
			}
			System.err.println("i_to_g.length: " + i_to_g.length);
			for (; ;) {
				boolean cut = true;
				String s = "";
				for (int i=0; i<phaseGenes.length; i++) {
					if (its[i].hasNext()) {
						s += old_i_to_g[its[i].next()] + "\t";
						cut = false;
					}
					else {
						s += "XXXX\t";
					}
				}
				if (cut) {
					break;
				}
				s = s.trim() + "\n";
				br.write(s);
			}
			br.close();
		}
		catch (IOException ie) {
			ie.printStackTrace();
		}
	}

	public void filterPaths() {
		System.err.println("l1, l2: " + l1 + ", " + l2);
		compute_initial_delta();
		/*
		for (int i=0; i<numGenes; ++i) {
			checkDelta(i);
		}
		*/
		int initialNumPaths = 0;
		double curObj = 0.;
		for (int i=0; i<numPaths; ++i) {
			initialNumPaths += b_p[i];
		}
		int initialNumGenes = 0;
		for (int i=0; i<numGenes; ++i) {
			initialNumGenes += b_g[i];
		}
		int initialNumTargets = 0;
		for (int i=0; i<index_to_targets.length; ++i) {
			initialNumTargets += b_g[index_to_targets[i]];
		}
		System.err.println("initialNumPaths: " + initialNumPaths);
		System.err.println("initialNumGenes: " + initialNumGenes);
		System.err.println("initialNumTargets: " + initialNumTargets);
		int numIterRun = 0;

		TIntLinkedList tabuList = new TIntLinkedList();
		Integer[] geneOrder = new Integer[numGenes];
		for (int i=0; i<geneOrder.length; ++i) {
			geneOrder[i] = i;
		}
		
		double bestSolutionObj = curObj;
		
		for (int i=0; i<b_g.length; ++i) {
			best_b_g[i] = b_g[i];
		}
		for (int i=0; i<b_p.length; ++i) {
			best_b_p[i] = b_p[i];
		}

		for (;numIterRun < N; ++numIterRun) {
			Arrays.sort(geneOrder, new GeneObjDeltaCompare(delta_g));

			/*
			double best_delta = -Double.MAX_VALUE;
			int g = -1;
			for (short i=0; i<numGenes; ++i) {
				if (delta_g[i] > best_delta) {
					best_delta = delta_g[i];
					g = i;
				}
			}
			System.err.println("best_delta:" + best_delta + ", " + 
					g + ", " + 
					initialNumGenes + ", " +
					initialNumTargets + ", " +
					b_g[g]);
			if (best_delta < 0.001 || b_g[g] == 0) { // TODO
				break;
			}
			*/

			int g = -1;
			for (int i=0; i<numGenes; ++i) {
				int curG = geneOrder[i];
				if (!tabuList.contains(curG)) {
					g = curG;
					break;
				}
			}
			System.err.println(numIterRun + ": bestG: " + g + ", " + delta_g[g] + ", " + b_g[g]);
			assert (g != -1);
			curObj += delta_g[g];
			
			boolean newBest = false;
			if (curObj > bestSolutionObj) {
				bestSolutionObj = curObj;
				newBest = true;
			}
			
			if (tabuList.size() >= L) {
				assert tabuList.size() == L : tabuList.size();
				tabuList.removeAt(0);
			}
			tabuList.add(g);
			
			double cache_delta2 = delta_g[g];
			add_target_penalty(-l2);
			double cache_delta = delta_g[g];
			if (b_g[g] == 1) {
				b_g[g] = 0;
				--initialNumGenes;
				delta_g[g] -= 2*l1;
				update_delta_on_deactivation(g);
			}
			else {
				++initialNumGenes;
				b_g[g] = 1;
				delta_g[g] += 2*l1;
				int numActivated = update_delta_on_activation(g);
				System.err.println("numActivated: " + numActivated);
				assert numActivated > 0 : numActivated;
			}
			if (newBest) {
				for (int i=0; i<b_g.length; ++i) {
					best_b_g[i] = b_g[i];
				}
				for (int i=0; i<b_p.length; ++i) {
					best_b_p[i] = b_p[i];
				}
			}
			
			assert Math.abs(delta_g[g]+cache_delta) < 0.1 : delta_g[g] + " + " +
			cache_delta + " == " + (delta_g[g]+cache_delta);
			double cache_delta3 = delta_g[g];
			initialNumTargets = add_target_penalty(l2);
			assert Math.abs(delta_g[g]+cache_delta2) < 0.1 : delta_g[g] + " + " +
			cache_delta2 + " == " + (delta_g[g]+cache_delta2) + ", " 
			+ cache_delta + " + " + cache_delta3;
		}
		System.err.println("numIterRun: " + numIterRun);

		for (int i=0; i<b_g.length; ++i) {
			b_g[i] = best_b_g[i];
		}
		for (int i=0; i<b_p.length; ++i) {
			b_p[i] = best_b_p[i];
		}
		
		int finalNumPaths = 0;
		for (int i=0; i<numPaths; ++i) {
			finalNumPaths += b_p[i];
		}
		int finalNumGenes = 0;
		for (int i=0; i<numGenes; ++i) {
			finalNumGenes += b_g[i];
		}
		int finalNumTargets = 0;
		for (int i=0; i<index_to_targets.length; ++i) {
			finalNumTargets += b_g[index_to_targets[i]];
		}
		System.err.println("finalNumPaths: " + finalNumPaths);
		System.err.println("finalNumGenes: " + finalNumGenes);
		System.err.println("finalNumTargets: " + finalNumTargets);
	}

	public void compute_initial_delta() {
		for (int i=0; i<numGenes; ++i) {
			delta_g[i] = l1;
			b_g[i] = 1;
		}
		for (int i=0; i<numPaths; ++i) {
			b_p[i] = 1;
			int target = P[i][P[i].length-1];
			int targetIndex = targets_to_index[target];
			f[targetIndex] += w_P[i];
			for (int j=0; j<P[i].length-1; ++j) {
				f_g_t[targetIndex][P[i][j]] -= w_P[i];
				delta_g[P[i][j]] -= w_P[i];
			}
			delta_g[target] -= w_P[i];
			f_g_t[targetIndex][target] -= w_P[i];
		}
		add_target_penalty(l2);
	}

	public int add_target_penalty(double penalty) {
		int numTargets = 0;
		for (int i=0; i<T.length; ++i) {
			//System.err.println(f[i] + ", " + T[i].length);
			if (f[i] >= 1) {
				++numTargets;
			}
			for (int j=0; j<T[i].length; ++j) {
				short g = T[i][j];
				if (f[i] >= 1) {
					if (f[i] + f_g_t[i][g] < 1) {
						delta_g[g] -= penalty;
					}
				}
				else {
					if (f[i] + f_g_t[i][g] >= 1) {
						delta_g[g] += penalty;
					}
				}
			}
		}
		return numTargets;
	}

	public int update_delta_on_activation(int g) {
		int num_activated = 0;
		for (int i=0; i<P_g[g].length; ++i) {
			int p = P_g[g][i];
			assert b_p[p] == 0;
			int r = P[p].length;
			short g3 = -1;
			for (int j=0; j<P[p].length; ++j) {
				short g2 = P[p][j];
				r -= b_g[g2];
				if (b_g[g2] == 0) {
					g3 = g2;
				}
			}
			int t = P[p][P[p].length-1];
			int tIndex = targets_to_index[t];
			
			// So the path needs to be activated now
			if (r == 0) {
				for (int j=0; j<P[p].length; ++j) {
					int g2 = P[p][j];
					delta_g[g2] -= w_P[p]; // because disabling g2 will disable the path now
					f_g_t[tIndex][g2] -= w_P[p];
				}
				f[tIndex] += w_P[p]; // flow increasing through target as path being activated
				delta_g[g] -= w_P[p]; // deactivating gene will lead to loss of path weight
				f_g_t[tIndex][g] -= w_P[p];
				b_p[p] = 1;
				++num_activated;
			}
			// The path is one gene flip away from activation
			else if (r == 1) {
				assert g3 != -1;
				assert g3 != g;
				delta_g[g3] += w_P[p];
				f_g_t[tIndex][g3] += w_P[p];
			}
		}
		return num_activated;
	}

	public void update_delta_on_deactivation(int g) {
		for (int i=0; i<P_g[g].length; ++i) {
			int p = P_g[g][i];
			b_p[p] = 0;
			int r = P[p].length;
			int inactive_gene = -1;
			for (int j=0; j<P[p].length; ++j) {
				short g2 = P[p][j];
				r -= b_g[g2];
				if (b_g[g2] == 0 && g != g2) {
					inactive_gene = g2;
				}
			}
			assert r >= 1;
			
			int t = P[p][P[p].length-1];
			int tIndex = targets_to_index[t];
			
			if (r == 2) {
				assert inactive_gene != -1;
				delta_g[inactive_gene] -= w_P[p];
				f_g_t[tIndex][inactive_gene] -= w_P[p];
			}
			else if (r == 1) {
				for (int j=0; j<P[p].length; ++j) {
					int g2 = P[p][j];
					delta_g[g2] += w_P[p];
					f_g_t[tIndex][g2] += w_P[p];
				}
				f[tIndex] -= w_P[p];
				delta_g[g] += w_P[p];
				f_g_t[tIndex][g] += w_P[p];
			}
		}
	}
	
	/*
	public int update_delta_on_activation(int g) {
		int num_activated = 0;
		for (int i=0; i<P_g[g].length; ++i) {
			int p = P_g[g][i];
			assert b_p[p] == 0;
			int r = P[p].length;
			for (int j=0; j<P[p].length; ++j) {
				short g2 = P[p][j];
				r -= b_g[g2];
			}
			int t = P[p][P[p].length-1];
			int tIndex = targets_to_index[t];

			// So the path needs to be activated now
			if (r == 0) {
				for (int j=0; j<P[p].length; ++j) {
					int g2 = P[p][j];
					delta_g[g2] -= w_P[p];
					f_g_t[tIndex][g2] -= w_P[p];
				}
				f[tIndex] += w_P[p];
				delta_g[g] -= w_P[p];
				f_g_t[tIndex][g] -= w_P[p];
				b_p[p] = 1;
				++num_activated;
			}
			// The path is one gene flip away from activation
			else if (r == 1) {
				delta_g[g] += w_P[p];
				f_g_t[tIndex][g] += w_P[p];
			}
		}
		return num_activated;
	}

	public void update_delta_on_deactivation(int g) {
		for (int i=0; i<P_g[g].length; ++i) {
			int p = P_g[g][i];
			int r = P[p].length;
			int inactive_gene = -1;
			for (int j=0; j<P[p].length; ++j) {
				short g2 = P[p][j];
				r -= b_g[g2];
				if (b_g[g2] == 0 && g != g2) {
					inactive_gene = g2;
				}
			}
			assert r >= 1;

			int t = P[p][P[p].length-1];
			int tIndex = targets_to_index[t];

			if (r == 2) {
				assert inactive_gene != -1;
				delta_g[inactive_gene] -= w_P[p];
				f_g_t[tIndex][inactive_gene] -= w_P[p];
			}
			else if (r == 1) {
				for (int j=0; j<P[p].length; ++j) {
					int g2 = P[p][j];
					delta_g[g2] += w_P[p];
					f_g_t[tIndex][g2] += w_P[p];
				}
				f[tIndex] -= w_P[p];
				delta_g[g] += w_P[p];
				f_g_t[tIndex][g] += w_P[p];
			}
		}
	}
	*/

	public void readData(HashMap<String, String> args) {
		g_to_i = new HashMap<String, Short>();

		BufferedReader br = null;
		try {
			HashSet<String> filteredGeneList = new HashSet<String>();
			String line;
			System.err.println("Read # of phases in");
			String phaseGenesFile = 
					argValue(args, "TPtargetsFile", "", true);
			br = new BufferedReader(new FileReader(phaseGenesFile));
			{
				line = br.readLine();
				String [] lineSplit = line.split("\\s+");
				numPhases = lineSplit.length;
				if (numPhases == 0) {
					System.err.println("ERROR: The file containing targets for each time point has a blank line at the top");
					System.exit(1);
				}

				phaseGenes = new TShortHashSet[numPhases];
				phaseTfs = new TShortHashSet[numPhases];
			}
			br.close();
			
			
			// Create initial filtered gene list from time series file
			System.err.println("Create initial filtered gene list from time series file");
			String timeseriesFilepath = 
					argValue(args, "timeseriesFilepath", "", true);
			br = new BufferedReader(new FileReader(timeseriesFilepath));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\t");
				filteredGeneList.add(lineSplit[0]);
				if (lineSplit[0].contains("GENE")) {
					continue;
				}
			}
			br.close();
			
			// Create gene name to index mapping from ppi while filtering
			// on filteredGeneList
			String ppiNetworkFilepath = 
					argValue(args, "ppiNetworkFilepath", "", true);
			br = new BufferedReader(new FileReader(ppiNetworkFilepath));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\\s+");
				String s = lineSplit[0].trim();
				String t = lineSplit[2].trim();

				if (!filteredGeneList.contains(s) 
						|| !filteredGeneList.contains(t)) 
				{
					continue;
				}

				if (!g_to_i.containsKey(s)) {
					assert g_to_i.size() <= Short.MAX_VALUE :
						Integer.toString(g_to_i.size());
					g_to_i.put(s, (short) g_to_i.size());
				}
				if (!g_to_i.containsKey(t)) {
					assert g_to_i.size() <= Short.MAX_VALUE
							: Integer.toString(g_to_i.size());
					g_to_i.put(t, (short) g_to_i.size());
				}
			}
			br.close();
			assert g_to_i.size() <= Short.MAX_VALUE
					: Integer.toString(g_to_i.size());

			numGenes = (short) g_to_i.size();
			System.err.println(numGenes + " genes in network");
			edgeScores = new TShortDoubleHashMap[numGenes];
			nl = new short[numGenes][];

			for (int i=0; i<numGenes; ++i) {
				edgeScores[i] = new TShortDoubleHashMap();
			}


			// Read in edges into edgeScores
			// So far every gene encountered will be in filtered gene list
			HashSet<String> ppi_incoming_edge_genes = new HashSet<String>();
			System.err.println("Reading in edges and their scores");
			int numEdges = 0;
			br = new BufferedReader(new FileReader(ppiNetworkFilepath));
			for (; (line = br.readLine()) != null; ++numEdges) {
				if (line == "") {
					continue;
				}
				String [] lineSplit = line.split("\\s+");
				assert lineSplit.length == 4 : line;
				String s = lineSplit[0].trim();
				String t = lineSplit[2].trim();

				if (!g_to_i.containsKey(s) ||
						!g_to_i.containsKey(t))
				{
					continue;
				}
				short tIndex = g_to_i.get(s);
				short sIndex = g_to_i.get(t);
				double score = Math.log(Double.parseDouble(lineSplit[3]));

				edgeScores[sIndex].put(tIndex, score);
				ppi_incoming_edge_genes.add(t);
				if (lineSplit[1].equals("(pp)")) {
					ppi_incoming_edge_genes.add(s);
					edgeScores[tIndex].put(sIndex, score);
				}
			}
			br.close();
			System.err.println(numEdges + " edges in network");


			// Create neighbor lists
			System.err.println("Create neighbor lists");
			class EdgeCompare implements Comparator<Short> {
				TShortDoubleHashMap edgeScores;
				public EdgeCompare(
						TShortDoubleHashMap edgeScores)
				{
					this.edgeScores = edgeScores;
				}

				public int compare(Short t1, Short t2) {
					double s1 = edgeScores.get(t1);
					double s2 = edgeScores.get(t2);
					// Descending order compare
					if (s1 < s2)
						return 1;
					else if (s1 > s2)
						return -1;
					else
						return 0;
				}
			}

			maxDegree = 0;
			for (short i=0; i<numGenes; ++i) {
				nl[i] = new short[edgeScores[i].size()];
				if (nl[i].length > maxDegree) {
					maxDegree = (short) nl[i].length;
				}
				Short [] temp = new Short[edgeScores[i].size()];
				int j=0;
				for (
						TShortDoubleIterator it = edgeScores[i].iterator();
						it.hasNext();
						++j)
				{
					it.advance();
					temp[j] = it.key();
				}
				assert j == temp.length;
				Arrays.sort(temp, new EdgeCompare(edgeScores[i]));
				if (DEBUG) {
					for (j=1; j<temp.length; ++j) {
						assert edgeScores[i].get(temp[j])
						<= edgeScores[i].get(temp[j-1]);
					}
				}
				for (j=0; j<temp.length; ++j) {
					nl[i][j] = temp[j];
				}
			}
			
			// Read Sources
			System.err.println("Read Sources");
			isSource = new boolean[numGenes];
			br = new BufferedReader(new FileReader(
					argValue(args, "sourcesFilepath", "", true)));
			int numSources = 0;
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\\s+");
				String source = lineSplit[0].trim();
				if (g_to_i.containsKey(source)) {
					++numSources;
					int sourceIndex = g_to_i.get(source);
					isSource[sourceIndex] = true;
				}
			}
			br.close();
			System.err.println(numSources + " sources in network");

			sources = new short[numSources];
			for (short i=0, j=0; i<numGenes; ++i) {
				if (isSource[i]) {
					sources[j++] = i;
				}
			}

			TShortHashSet tfsWithIncomingInteractionInfo = new TShortHashSet();
			TShortHashSet genesWithInteractionInfo = new TShortHashSet();
			String tfGeneFile = 
					argValue(args, "tfGeneFile", "", true);
			br = new BufferedReader(new FileReader(tfGeneFile));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\\s+");
				if (g_to_i.containsKey(lineSplit[0])
						&& g_to_i.containsKey(lineSplit[1])) 
				{
					if (ppi_incoming_edge_genes.contains(lineSplit[0])
							|| isSource[g_to_i.get(lineSplit[0])])
					{
						genesWithInteractionInfo.add(g_to_i.get(lineSplit[1]));
						tfsWithIncomingInteractionInfo.add(g_to_i.get(lineSplit[0]));
					}
				}
			}
			br.close();

			// Read phase genes
			for (int i=0; i<numPhases; ++i) {
				System.err.println("Reading genes for phase " + i);
				br = new BufferedReader(new FileReader(phaseGenesFile));
				phaseGenes[i] = new TShortHashSet();
				for (; (line = br.readLine()) != null;) {
					String [] lineSplit = line.split("\t");
					assert phaseGenes[i].size() <= numPhaseGenes 
							: phaseGenes[i].size();
					if (lineSplit[i] == "XXXX") {
						continue;
					}
					if (phaseGenes[i].size() == numPhaseGenes) {
						break;
					}
					if (g_to_i.containsKey(lineSplit[i])) {
						short gene = g_to_i.get(lineSplit[i]);
						if (!genesWithInteractionInfo.contains(gene)) {
							continue;
						}
						boolean useGene = true;
						for (int j=0; j<i; ++j) {
							if (phaseGenes[j].contains(gene)) {
								useGene = false;
								break;
							}
						}
						if (useGene) {
							phaseGenes[i].add(gene);
						}
					}
				}
				System.err.println(phaseGenes[i].size() 
						+ " phase genes for phase " + i);
				br.close();
			}
			
			String phaseTfFilepath = argValue(args, "phaseTfFilepath", "", false);
			if (phaseTfFilepath != "") {
				// Read phase Tfs
				for (int i=0; i<numPhases; ++i) {
					System.err.println("Reading Tfs for phase " + i);
					br = new BufferedReader(new FileReader(phaseTfFilepath));
					phaseTfs[i] = new TShortHashSet();
					for (; (line = br.readLine()) != null;) {
						String [] lineSplit = line.split("\t");
						if (lineSplit[i].startsWith("XXXX")) {
							continue;
						}
						if (g_to_i.containsKey(lineSplit[i])) {
							short tf = g_to_i.get(lineSplit[i]);
							if (!tfsWithIncomingInteractionInfo.contains(tf)) {
								continue;
							}
							phaseTfs[i].add(tf);
						}
					}
					System.err.println(phaseTfs[i].size() 
							+ " phase Tfs for phase " + i);
					br.close();
				}
			}
			
			/*
			// Compute fold change
			foldChange = new double[numGenes][];
			br = new BufferedReader(new FileReader(timeseriesFilepath));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\t");
				if (!g_to_i.containsKey(lineSplit[0])) {
					continue;
				}
				short gene = g_to_i.get(lineSplit[0]);
				double [] expr;
				int numTp = lineSplit.length-1;
				double [] change = new double[numPhases];
				assert numTp == numPhases;
				for (int i=1; i<lineSplit.length; ++i) {
					change[i-1] = Double.parseDouble(lineSplit[i]);
				}
				foldChange[gene] = change;
			}
			br.close();
			*/

			// Read TF-gene interactions
			System.err.println("Reading TFs");
			tfToGene = new TShortDoubleHashMap[numGenes][numPhases];
			isTf = new boolean[numGenes];
			isMirna = new boolean[numGenes];
			isRnaHit = new boolean[numGenes];
			String rnaHitsFile = argValue(args, "rnaHitsFile", "", false);
			for (int i=0; i<numGenes; ++i) {
				isTf[i] = false;
				isMirna[i] = false;
				isRnaHit[i] = false;
				for (int j=0; j<numPhases; ++j) {
					tfToGene[i][j] = new TShortDoubleHashMap();
				}
			}
			br = new BufferedReader(new FileReader(tfGeneFile));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\\s+");
				if (g_to_i.containsKey(lineSplit[0]) 
						&& tfsWithIncomingInteractionInfo.contains(
								g_to_i.get(lineSplit[0]))) 
				{
					short tfIndex = g_to_i.get(lineSplit[0]);
					isMirna[tfIndex] = lineSplit[0].startsWith(mirnaPrefix);
					if (g_to_i.containsKey(lineSplit[1])) 
					{
						short geneIndex = g_to_i.get(lineSplit[1]);
						double score = Math.log(Double.parseDouble(lineSplit[2]));
						for (int j=0; j<numPhases; ++j) {
							if (phaseGenes[j].contains(geneIndex)) {
								isTf[tfIndex] = true;
								tfToGene[tfIndex][j].put(geneIndex, score);
							}
						}
					}
				}
			}
			br.close();
			br = new BufferedReader(new FileReader(tfGeneFile));
			for (; (line = br.readLine()) != null;) {
				String [] lineSplit = line.split("\\s+");
				if (g_to_i.containsKey(lineSplit[0])) {
					isRnaHit[g_to_i.get(lineSplit[0])] = true;
				}
			}
			br.close();
			int numTfs = 0;
			for (int i=0; i<numGenes; ++i) {
				numTfs += isTf[i] ? 1 : 0;
			}
			System.err.println(numTfs + " tfs in network");
		}
		catch(IOException ie) {
			ie.printStackTrace();
		}

		pathBufferSize = numGenes > pathBufferSize ? numGenes : pathBufferSize;
		i_to_g = new String[g_to_i.size()];

		// Map gene index to name
		for (Map.Entry<String, Short> entry : g_to_i.entrySet()) {
			i_to_g[entry.getValue()] = entry.getKey();
		}
	}
}
