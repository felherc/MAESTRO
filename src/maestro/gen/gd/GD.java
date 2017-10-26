/**
 * Copyright 2017 Felipe Hernández
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the 
 * License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 * express or implied. See the License for the specific language governing permissions and
 * limitations under the License.
 */

package maestro.gen.gd;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeSet;

import probDist.Normal;
import probDist.Uniform;
import utilities.Utilities;
import utilities.stat.ContSeries;
import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.Generator;
import maestro.gen.gd.IndexDoubleSorter;
import maestro.gen.gd.NormSolution;
import maestro.gen.gd.Vicinity;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;

/**
 * This generator allows to create new solutions to a multi-objective optimization problem with 
 * discrete and continuous decision variables using a population-based Gradient Descent (GD) 
 * algorithm. Solutions are selected from the population and a vicinity is created around each of 
 * them with the other solutions until there are enough to compute empirical steepest gradients for 
 * all the objectives. New solutions then are sampled in the direction of a random linear 
 * combination of the gradient vectors.
 * @author Felipe Hernández
 */
public class GD implements Generator
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator method
	 */
	public final static String ID = "Gradient Descent";
	
	/**
	 * The short version of the identifier of the generator method
	 */
	public final static String SHORT_ID = "GD";
	
	/**
	 * Report name for {@link #greed}
	 */
	public static final String PARAM_GREED = "Greed";
	
	/**
	 * Report name for {@link #solVicinities}
	 */
	public static final String PARAM_SOL_VICINITIES = "Solutions in vicinities";
	
	/**
	 * Report name for {@link #powerVicinityWeight}
	 */
	public static final String PARAM_POWER_VICINITY_WEIGHT = "Vicinity weight power";
	
	/**
	 * Report name for {@link #stepSize}
	 */
	public static final String PARAM_STEP_SIZE = "Step size";
	
	/**
	 * Report name for {@link #amplitude}
	 */
	public static final String PARAM_AMPLITUDE = "Amplitude";
	
	/**
	 * Report name for {@link #uniformProb}
	 */
	public static final String PARAM_UNIFORM_PROB = "Uniform probability";
	
	/**
	 * Minimum value for the selection weight assignment variable <i>q</i>, corresponding to a 
	 * maximum weight for the solutions of rank one. Weights are assigned according to the 
	 * probability density of a Normal distribution with mean 1.0 and standard deviation <i>qk</i>, 
	 * where <i>k</i> is the number of solutions in the population. The argument for the 
	 * distribution is the rank of the solution.
	 */
	public final static double MIN_Q = 0.1;
	
	/**
	 * Maximum value for the selection weight assignment variable <i>q</i>, corresponding to a 
	 * uniform probability for all the solutions. Weights are assigned according to the probability 
	 * density of a Normal distribution with mean 1.0 and standard deviation <i>qk</i>, where 
	 * <i>k</i> is the number of solutions in the population. The argument for the distribution is 
	 * the rank of the solution.
	 */
	public final static double MAX_Q = 10.0;
	
	/**
	 * Power of the function to obtain the value of <i>q</i> for weight assignment from different
	 * greed values. A higher power increases the greed range in which the solutions in the first
	 * front are preferred; a lower power increases the range in which all solutions are equally
	 * likely to be selected.
	 */
	public final static double GREED_TO_Q_POW = 5.0;

	/**
	 * The maximum number of attempts to generate the each candidate base for the vicinities in 
	 * terms of the ratio of the number of required bases. Attempts fail if the sampled base was
	 * sampled before.
	 */
	private final static double MAX_SAMPLE_ATTEMPTS = 3.0;
	
	/**
	 * Default value for {@link #greed}
	 */
	public final static double DEF_GREED = 0.75;
	
	/**
	 * Default value for {@link #solVicinities}
	 */
	public final static double DEF_SOL_VICINITIES = 1.0;
	
	/**
	 * Default value for {@link #powerVicinityWeight}
	 */
	public final static double DEF_POWER_VICINITY_WEIGHT = 0.0;
	
	/**
	 * Default value for {@link #stepSize}
	 */
	public final static double DEF_STEP_SIZE = 0.0001;
	
	/**
	 * Default value for {@link #amplitude}
	 */
	public final static double DEF_AMPLITUDE = 1.5;
	
	/**
	 * Default value for {@link #uniformProb}
	 */
	public final static double DEF_UNIFORM_PROB = 0.0;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The list of the discrete decision variables of the problem
	 */
	private ArrayList<DiscVar> discVars;
	
	/**
	 * The list of the continuous decision variables of the problem
	 */
	private ArrayList<ContVar> contVars;
	
	/**
	 * A value between -1.0 and 1.0 that represents the expected mean "quality" of the solutions
	 * to be selected from the current population to serve as the base of the vicinities. 1.0 if 
	 * the selected solutions should only be among the best according to the objectives. -1.0 if 
	 * they should only be among the worst. 0.0 if there is no preference and the solutions have 
	 * uniform chances to be selected.
	 */
	private double greed;
	
	/**
	 * The percentage of the number of solutions in the parent population to be used for the 
	 * creation of vicinities. For example, if there are nine solutions and there are two decision
	 * variables, a value of 2.0 would mean that six vicinities would be created, each of which of 
	 * three solutions, to account for two times the number of solutions in the population.
	 */
	private double solVicinities;
	
	/**
	 * The power to apply to the magnitude of the gradient sum of each vicinity to be used as the
	 * weight to compute the probability of using a specific vicinity to generate new solutions.
	 * That is, the probability of generating solutions from a vicinity is proportional to the 
	 * magnitude of the sum of all its gradients to the power specified.
	 */
	private double powerVicinityWeight;
	
	/**
	 * Also called learning rate, it is a factor to multiply the gradient vector of a vicinity to 
	 * determine the distance from the vicinity's base to create new solutions. The resulting point
	 * is used as the center of a Gaussian kernel from which new solutions are sampled. Small 
	 * values for the step size would account for higher accuracy in finding optima, but would also 
	 * require more iterations to reach them. Large values can accelerate the search but can also 
	 * make convergence harder as the large steps would tend to overshoot the optima once it is 
	 * close.
	 */
	private double stepSize;
	
	/**
	 * This factor defines the breadth of the kernel to sample new solutions from. The factor
	 * multiplies each of the values in the diagonal of the covariance matrix of the Gaussian 
	 * kernel, which are assigned the variance of the vicinity on each dimension by default. Large 
	 * values for the amplitude would make the sampled solutions more likely to drift away from the 
	 * direction of the gradients.
	 */
	private double amplitude;
	
	/**
	 * The probability of each value of a generated solution to be re-sampled from a uniform 
	 * distribution spanning the entire valid range of the variable. Analogous to the mutation
	 * probability in Genetic Algorithms.
	 */
	private double uniformProb;
		
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new instance of the Gradient Descent generator with the default parameters
	 */
	public GD()
	{
		discVars			= new ArrayList<>();
		contVars			= new ArrayList<>();
		greed				= DEF_GREED;
		solVicinities		= DEF_SOL_VICINITIES;
		powerVicinityWeight	= DEF_POWER_VICINITY_WEIGHT;
		stepSize			= DEF_STEP_SIZE;
		amplitude			= DEF_AMPLITUDE;
		uniformProb			= DEF_UNIFORM_PROB;
	}
	
	/**
	 * Creates a new instance of the Gradient Descent generator
	 * @param greed					{@link #greed}
	 * @param solVicinities			{@link #solVicinities}
	 * @param powerVicinityWeight	{@link #powerVicinityWeight}
	 * @param stepSize				{@link #stepSize}
	 * @param amplitude				{@link #amplitude}
	 * @param uniformProb			{@link #uniformProb}
	 */
	public GD(double greed, double solVicinities, double powerVicinityWeight, double stepSize, 
				double amplitude, double uniformProb) 
	{
		discVars					= new ArrayList<>();
		contVars					= new ArrayList<>();
		this.greed					= greed;
		this.solVicinities			= solVicinities;
		this.powerVicinityWeight	= powerVicinityWeight;
		this.stepSize				= stepSize;
		this.amplitude				= amplitude;
		this.uniformProb			= uniformProb;
	}

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	@Override
	public String getId() 
	{
		return ID;
	}

	@Override
	public String getShortId() 
	{
		return SHORT_ID;
	}

	/**
	 * @return {@link #greed}
	 */
	public double getGreed() 
	{
		return greed;
	}

	/**
	 * @param greed {@link #greed}
	 */
	public void setGreed(double greed) 
	{
		this.greed = greed;
	}

	/**
	 * @return {@link #solVicinities}
	 */
	public double getSolVicinities() 
	{
		return solVicinities;
	}

	/**
	 * @param solVicinities {@link #solVicinities}
	 */
	public void setSolVicinities(double solVicinities) 
	{
		this.solVicinities = solVicinities;
	}

	/**
	 * @return {@link #powerVicinityWeight}
	 */
	public double getPowerVicinityWeight() 
	{
		return powerVicinityWeight;
	}

	/**
	 * @param powerVicinityWeight {@link #powerVicinityWeight}
	 */
	public void setPowerVicinityWeight(double powerVicinityWeight) 
	{
		this.powerVicinityWeight = powerVicinityWeight;
	}

	/**
	 * @return {@link #stepSize}
	 */
	public double getStepSize() 
	{
		return stepSize;
	}

	/**
	 * @param stepSize {@link #stepSize}
	 */
	public void setStepSize(double stepSize) 
	{
		this.stepSize = stepSize;
	}

	/**
	 * @return {@link #amplitude}
	 */
	public double getAmplitude() 
	{
		return amplitude;
	}

	/**
	 * @param amplitude {@link #amplitude}
	 */
	public void setAmplitude(double amplitude) 
	{
		this.amplitude = amplitude;
	}
	
	/**
	 * @return {@link #uniformProb}
	 */
	public double getUniformProb() 
	{
		return uniformProb;
	}

	/**
	 * @param uniformProb {@link #uniformProb}
	 */
	public void setUniformProb(double uniformProb) 
	{
		this.uniformProb = uniformProb < 0.0 ? 0.0 : uniformProb;
	}

	@Override
	public String getParamSummary() 
	{
		String line	= "";
		line		+=			PARAM_GREED					+ " = " + greed;
		line		+= "; " +	PARAM_SOL_VICINITIES		+ " = " + solVicinities;
		line		+= "; " +	PARAM_POWER_VICINITY_WEIGHT	+ " = " + powerVicinityWeight;
		line		+= "; " +	PARAM_STEP_SIZE				+ " = " + stepSize;
		line		+= "; " +	PARAM_AMPLITUDE				+ " = " + amplitude;
		line		+= "; " +	PARAM_UNIFORM_PROB			+ " = " + uniformProb;
		return line;
	}

	@Override
	public void addDiscVariable(DiscVar variable) 
	{
		discVars.add(variable);
	}

	@Override
	public void addContVariable(ContVar variable) 
	{
		contVars.add(variable);
	}

	@Override
	public void clearVariables() 
	{
		discVars.clear();
		contVars.clear();
	}

	@Override
	public ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population,
														int number)
	{
		// Compute statistics for the scalar decision variables
		int discVarCount					= discVars.size();
		int contVarCount					= contVars.size();
		ArrayList<Integer> discScalarIDs	= new ArrayList<>();
		ArrayList<ContSeries> scalarStats	= new ArrayList<>();
		for (int j = 0; j < contVarCount; j++)
			scalarStats.add(new ContSeries(false));
		for (int j = 0; j < discVarCount; j++)
		{
			DiscVar discVar 				= discVars.get(j);
			if (discVar.isScalar())
			{
				discScalarIDs.add(j);
				scalarStats.add(new ContSeries(false));
			}
		}
		for (Solution solution : population)
		{
			ArrayList<Double> contValues	= solution.getContValues();
			for (int j = 0; j < contVars.size(); j++)
				scalarStats.get(j).addValue(contValues.get(j));
			ArrayList<Integer> discValues	= solution.getDiscValues();
			for (Integer j : discScalarIDs)
				scalarStats.get(j + contVarCount).addValue(discValues.get(j));
		}
		
		// Verify if there are zero-variance dimensions and if there are enough solutions to
		// compute the gradients
		HashSet<Integer> zeroVar			= new HashSet<>();
		for (int v = 0; v < scalarStats.size(); v++)
			if (scalarStats.get(v).getVar() == 0.0)
				zeroVar.add(v);
		int validScalarCount				= scalarStats.size() - zeroVar.size();
		if (validScalarCount > population.size() - 1)
			return new ArrayList<>();
		
		// Create normalized solutions
		ArrayList<NormSolution> normSolList	= createNormalizedSolutions(contVarCount, 
												discScalarIDs, scalarStats, zeroVar, population);
		
		// Select bases for vicinities
		TreeSet<Integer> baseIDs			= getVicinityBases(population, validScalarCount);
		
		// Compute distance matrix
		ArrayList<TreeSet<IndexDoubleSorter>> sortings = computeDistances(normSolList, baseIDs);
		
		// Create vicinities
		ArrayList<Vicinity> vicinities		= new ArrayList<>();
		Iterator<Integer> iterIDs			= baseIDs.iterator();
		for (int i = 0; i < baseIDs.size(); i++)
		{
			int baseIndex					= iterIDs.next();
			NormSolution base				= normSolList.get(baseIndex);
			base.setBaseIndex(baseIndex);
			Vicinity vicinity				= new Vicinity(base);
			Iterator<IndexDoubleSorter> iter = sortings.get(i).iterator();
			iter.next();	// Discard first (the base itself)
			while (vicinity.size() < validScalarCount + 1 && iter.hasNext())
			{
				int index					= iter.next().getIndex();
				NormSolution candidate		= normSolList.get(index);
				vicinity.offerSolution(candidate);
			}
			if (vicinity.size() >= validScalarCount + 1)
				vicinities.add(vicinity);
		}
		
		// Determine the number of solutions to be generated by each vicinity
		ContSeries selector					= new ContSeries(true);
		int[] toGenerate					= new int[vicinities.size()];
		for (int v = 0; v < vicinities.size(); v++)
		{
			Vicinity vicinity				= vicinities.get(v);
			double gradMagnitude			= vicinity.getGradMagnitude();
			if (!Double.isNaN(gradMagnitude))
				selector.addValue(v, Math.pow(gradMagnitude, powerVicinityWeight));
			toGenerate[v]					= 0;
		}
		ArrayList<Double> vicinityIndexes	= selector.sampleMultiple(number, true, true);
		if (vicinityIndexes == null)
			return new ArrayList<>();
		for (Double index : vicinityIndexes)
			toGenerate[(int)(double)index]++;
			
		// Generate solutions
		ArrayList<Vicinity> activeVicin		= new ArrayList<>(vicinities);
		ArrayList<NormSolution> solutions	= new ArrayList<>();
		int index 							= 0;
		int invalidCount					= 0;
		while (solutions.size() < number && invalidCount < activeVicin.size())
		{
			Vicinity vicinity				= activeVicin.get(index);
			if (vicinity != null)
			{
				int toGen = Math.min(toGenerate[index], number - solutions.size());
				if (toGen > 0)
				{
					ArrayList<NormSolution> cand = vicinity.randomGenerate(toGen, stepSize, 
																			amplitude);
					if (cand != null)
						solutions.addAll(cand);
					else
					{
						activeVicin.set(index, null);
						invalidCount++;
					}
				}
			}
			index++;
			if (index >= activeVicin.size())
				index = 0;
		}
		
		// Un-normalize generated solutions
		ArrayList<SolutionRoot> roots		= unnormalizeSolutions(discVarCount, contVarCount, 
									discScalarIDs, scalarStats, zeroVar, population, solutions);
		
		// Apply uniform re-sampling
		if (uniformProb > 0.0)
			for (SolutionRoot root : roots)
				for (int i = 0; i < discVarCount + contVarCount; i++)
					if (Math.random() < uniformProb)
					{
						if (i < discVarCount)
						{
							DiscVar var	= discVars.get(i);
							int value	= Utilities.uniformRandomSelect(var.getCount()) 
											+ var.getMin();
							root.getDiscValues().set(i, value);
						}
						else
						{
							index		= i - discVarCount;
							ContVar var	= contVars.get(index);
							double val	= Uniform.sample(var.getMin(), var.getMax());
							root.getContValues().set(index, val);
						}
					}
		
		// Return solution roots
		return roots;
	}

	/**
	 * Computes the distance between the bases and all the solutions and returns the sets of indices
	 * of the solutions sorted by distance
	 * @param normSolList The list of normalized solutions
	 * @param baseIDs The indexes of the base solutions
	 * @return The sets of indices of the solutions sorted by distance
	 */
	private ArrayList<TreeSet<IndexDoubleSorter>> computeDistances(
			ArrayList<NormSolution> normSolList, TreeSet<Integer> baseIDs) 
	{
		ArrayList<TreeSet<IndexDoubleSorter>> sortings = new ArrayList<>();
		for (Integer baseID : baseIDs)
		{
			NormSolution base				= normSolList.get(baseID);
			TreeSet<IndexDoubleSorter> sort = new TreeSet<>();
			for (int s = 0; s < normSolList.size(); s++)
			{
				double distance				= Double.NaN;
				if (s == baseID)
					distance				= 0.0;
				else
				{
					NormSolution other		= normSolList.get(s);
					distance				= base.computeDistance(other);
				}
				sort.add(new IndexDoubleSorter(s, distance));
			}
			sortings.add(sort);
		}
		return sortings;
	}

	/**
	 * Normalizes the values for the decision variables and the fitness for the list of solutions
	 * @param contVarCount		The number of continuous decision variables
	 * @param discScalarIDs		The indexes of the scalar discrete decision variables
	 * @param scalarStats 		The list of statistics for the scalar decision variables
	 * @param zeroVar			A set containing the indices of the scalar decision variables whose
	 * 							variance is zero. These variables are not included in the 
	 * 							normalized solutions.
	 * @param solList			The list of unnormalized solutions
	 * @return The list of solutions with normalized values
	 */
	private ArrayList<NormSolution> createNormalizedSolutions(int contVarCount,
			ArrayList<Integer> discScalarIDs, ArrayList<ContSeries> scalarStats,
			HashSet<Integer> zeroVar, ArrayList<Solution> solList)
	{
		ArrayList<NormSolution> normSolList = new ArrayList<>();
		int scalarCount						= scalarStats.size();
		for (int s = 0; s < solList.size(); s++)
		{
			// Normalize values for decision variables
			Solution solution 				= solList.get(s);
			double[] normValues				= new double[scalarCount - zeroVar.size()];
			ArrayList<Double> contValues	= solution.getContValues();
			ArrayList<Integer> discValues	= solution.getDiscValues();
			int targetIndex					= 0;
			for (int j = 0; j < scalarCount; j++)
			{
				if (!zeroVar.contains(new Integer(j)))
				{
					if (j < contVarCount)
					{
						double original			= contValues.get(j);
						ContSeries stats		= scalarStats.get(j);
						normValues[targetIndex]	= (original - stats.getMean())/stats.getStDev();
					}
					else
					{
						int index				= discScalarIDs.get(j - contVarCount);
						double original			= discValues.get(index);
						DiscVar var				= discVars.get(index);
						Uniform dist			= new Uniform(var.getMin(), var.getMax());
						normValues[targetIndex]	= (original - dist.getMean())/dist.getStDev();
					}
					targetIndex++;
				}
			}
			
			// Normalize fitness value
			Uniform dist			= new Uniform(0, solList.size() - 1);
			double normFitness		= (s - dist.getMean())/dist.getStDev();

			normSolList.add(new NormSolution(normValues, normFitness));
		}
		return normSolList;
	}

	/**
	 * Selects <i>n</i> bases of the vicinities. <i>n</i> depends on the number of solutions,
	 * the number of decision variables, and {@link #solVicinities}. The bases are selected based
	 * on their rank and on {@link #greed}.
	 * @param solList The list with the solutions in the current population
	 * @param scalarCount The number of scalar decision variables
	 * @return The set of 
	 */
	private TreeSet<Integer> getVicinityBases(ArrayList<Solution> solList, int scalarCount) 
	{		
		// Determine the number of vicinities to create
		TreeSet<Integer> baseIDs	= new TreeSet<>();
		int vicinityCount			= (int)Math.max(1.0, 
										solVicinities*solList.size()/(scalarCount + 1));
		if (vicinityCount >= solList.size())
		{
			for (int i = 0; i < solList.size(); i++)
				baseIDs.add(i);
			return baseIDs;
		}
		
		// Obtain q
		greed						= greed >  1.0 ?  1.0 : greed;
		greed						= greed < -1.0 ? -1.0 : greed;		
		double temp					= Math.pow(1 - Math.abs(greed), GREED_TO_Q_POW);
		double q					= MIN_Q + (MAX_Q - MIN_Q)*temp;
		
		// Assign weights
		ContSeries selector			= new ContSeries(true);
		double stDev				= q*solList.size();
		for (int i = 0; i < solList.size(); i++)
		{
			double weight			= Normal.computepdf(1.0, stDev, i);
			selector.addValue(i, weight);
		}
		
		// Select bases for vicinities
		int maxAttempts				= (int)(vicinityCount*MAX_SAMPLE_ATTEMPTS);
		int attempts				= 0;
		while (baseIDs.size() < vicinityCount && attempts <= maxAttempts)
		{
			baseIDs.add((int)selector.sample());
			attempts++;
		}
		return baseIDs;
	}

	/**
	 * Un-normalizes the list of solutions to the original space
	 * @param discVarCount	The number of discrete decision variables
	 * @param contVarCount	The number of continuous decision variables
	 * @param discScalarIDs	The indexes of the scalar discrete decision variables
	 * @param scalarStats	The list of statistics for the scalar decision variables
	 * @param zeroVar		A set containing the indices of the scalar decision variables whose
	 * 						variance is zero. These variables are not included in the 
	 * 						normalized solutions.
	 * @param solList		The list with the solutions in the current population
	 * @param solutions		The list of normalized solutions
	 * @return The list of un-normalized solutions
	 */
	private ArrayList<SolutionRoot> unnormalizeSolutions(int discVarCount, int contVarCount,
			ArrayList<Integer> discScalarIDs, ArrayList<ContSeries> scalarStats,
			HashSet<Integer> zeroVar, ArrayList<Solution> solList, 
			ArrayList<NormSolution> solutions)
	{
		ArrayList<SolutionRoot> roots		= new ArrayList<>();
		for (NormSolution normSolution : solutions)
		{
			ArrayList<Double> contValues	= contVarCount > 0 ? new ArrayList<>() : null;
			ArrayList<Integer> discValues	= discVarCount > 0 ? new ArrayList<>() : null;
			Solution base					= solList.get(normSolution.getBaseIndex());
			int normIndex					= 0;
			for (int v = 0; v < scalarStats.size(); v++)
			{
				if (!zeroVar.contains(v))
				{
					double normalized		= normSolution.getValues()[normIndex];
					ContSeries stats		= scalarStats.get(v);
					double unnormalized		= normalized*stats.getStDev() + stats.getMean();
					if (v < contVarCount)
						contValues.add(unnormalized);
					else
					{
						int index			= discScalarIDs.get(v - contVarCount);
						discValues.set(index, discVars.get(index).validate((int)unnormalized));
					}
					normIndex++;
				}
				else
				{
					if (v < contVarCount)
						contValues.add(base.getContValues().get(v));
					else
					{
						int index			= discScalarIDs.get(v - contVarCount);
						discValues.set(index, base.getDiscValues().get(v));
					}
				}
			}
			int discScalarIndex				= 0;
			for (int v = 0; v < discVarCount; v++)
			{
				if (v != discScalarIDs.get(discScalarIndex))
					discValues.set(v, base.getDiscValues().get(v));			
				else
					discScalarIndex++;
			}
			roots.add(new SolutionRoot(discValues, contValues));
		}
		return roots;
	}

}