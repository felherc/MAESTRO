package maestro.gen.aco;

import java.util.ArrayList;
import java.util.Collections;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.Generator;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;
import probDist.ContProbDist;
import probDist.DiscProbDist;
import probDist.KernelDensity;
import probDist.Normal;
import probDist.Uniform;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.MultiVarNormal;
import probDist.multiVar.tools.Sample;
import utilities.MatUtil;
import utilities.Utilities;
import utilities.stat.ContSeries;

/**
 * This generator allows to create new solutions to a problem with discrete and continuous 
 * variables using a hybrid approach between a Metropolis random-walk adaptive sampling
 * algorithm, and an Ant Colony Optimization (ACO) algorithm. The ACO approach creates new 
 * solutions by sampling values independently from marginalized distributions of the current 
 * population for each variable. The Metropolis approach samples the values from a joint
 * distribution instead. The hybrid approach first generates a subset of values using the ACO 
 * approach, and then generates the rest using the Metropolis approach on a distribution 
 * conditioned on the ACO values.
 * @author Felipe Hernández
 */
public class MetroACO implements Generator 
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator method
	 */
	public final static String ID = "Metropolis - Ant Colony Optimization";
	
	/**
	 * The short version of the identifier of the generator method
	 */
	public final static String SHORT_ID = "MetroACO";
	
	/**
	 * Report name for {@link #greed}
	 */
	public static final String PARAM_GREED = "Greed";
	
	/**
	 * Report name for {@link #independentRatio}
	 */
	public static final String PARAM_INDEPENDENT_RATIO = "Independent ratio";
	
	/**
	 * Report name for {@link #uniformProb}
	 */
	public static final String PARAM_UNIFORM_PROB = "Uniform probability";
	
	/**
	 * Report name for {@link #groupNumber}
	 */
	public static final String PARAM_GROUP_NUMBER = "Group number";
	
	/**
	 * Report name for {@link #groupPercent}
	 */
	public static final String PARAM_GROUP_PERCENT = "Group percentage";
	
	/**
	 * Report name for {@link #bandwidthMethod}
	 */
	public static final String PARAM_BANDWIDTH_METHOD = "Bandwidth method";
	
	/**
	 * Report name for {@link #bandwidthMult}
	 */
	public static final String PARAM_BANDWIDTH_MULT = "Bandwidth multiplier";
	
	/**
	 * Report name for {@link #BANDWIDTH_SILVERMAN_UNIVAR}
	 */
	public static final String PARAM_B_W_SILVERMAN_UNIVAR = "Univariate Silverman";
	
	/**
	 * Report name for {@link #BANDWIDTH_SCOTT}
	 */
	public static final String PARAM_B_W_SCOTT = "Scott";
	
	/**
	 * Report name for {@link #BANDWIDTH_SILVERMAN}
	 */
	public static final String PARAM_B_W_SILVERMAN = "Silverman";
	
	/**
	 * {@link #bandwidthMethod} option: Silverman's rule for univariate distributions: the rule is
	 * applied on each dimension individually. The coefficient for the bandwidth (<i>h</i>) matrix 
	 * is computed as: <i>h = [(4*sig^5)/(3*n)]^(2/5)</i>, where <i>sig</i> is the standard
	 * deviation of the samples for that dimension, and <i>n</i> is the number of samples.
	 */
	public final static int BANDWIDTH_SILVERMAN_UNIVAR = 0;
	
	/**
	 * {@link #bandwidthMethod} option: Scott's rule: Estimates the bandwidth of each dimension 
	 * of the kernel density distribution of the root solutions using the optimal criteria for an 
	 * assumed underlying Gaussian distributions as described by Scott's rule.
	 */
	public final static int BANDWIDTH_SCOTT = 1;
	
	/**
	 * {@link #bandwidthMethod} option: Silverman's rule: Estimates the bandwidth of each dimension 
	 * of the kernel density distribution of the root solutions using the optimal criteria for an 
	 * assumed underlying Gaussian distributions as described in the reference below. The criteria 
	 * is termed "normal distribution approximation" or "Silverman's rule of thumb".<br><br>
	 * Silverman, B.W., 1998, <i>Density Estimation for Statistics and Data Analysis</i>, London: 
	 * Chapman & Hall/CRC, p.48, ISBN 0-412-24620-1.
	 */
	public final static int BANDWIDTH_SILVERMAN = 2;
	
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
	 * Default value for {@link #greed}
	 */
	public final static double DEF_GREED = 0.9;
	
	/**
	 * Default value for {@link #independentRatio}
	 */
	public final static ContProbDist DEF_INDEPENDENT_RATIO = new Uniform(0.0, 0.0);
	
	/**
	 * Default value for {@link #uniformProb}
	 */
	public final static double DEF_UNIFORM_PROB = 0.025;
	
	/**
	 * Default value for {@link #groupNumber}
	 */
	public final static int DEF_GROUP_NUMBER = -1;
	
	/**
	 * Default value for {@link #groupPercent}
	 */
	public final static double DEF_GROUP_PERCENT = 1.5;
	
	/**
	 * Default value for {@link #bandwidthMethod}
	 */
	public final static int DEF_BANDWIDTH_METHOD = BANDWIDTH_SILVERMAN;
	
	/**
	 * Default value for {@link #bandwidthMult}
	 */
	public final static double DEF_BANDWIDTH_MULT = 1.0;
	
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
	 * to be selected from the current population. 1.0 if the selected solutions should only be 
	 * among the best according to the objectives. -1.0 if they should only be among the worst. 
	 * 0.0 if there is no preference and the solutions have uniform chances to be selected.
	 */
	private double greed;
	
	/**
	 * Probability distribution for selecting the ratio of variables (<i>p</i>) that are to be 
	 * sampled independently, in an ACO fashion, for the creation of new solutions. The values are 
	 * sampled from a marginalized probability distribution for a number of the variables that 
	 * are selected randomly. The ratio of the rest of the values of the new solution 
	 * <i>(1 - p)</i> are sampled from a distribution conditioned on the independently sampled 
	 * values. For each new solution to be created, <i>p</i> is sampled from this distribution. If 
	 * the value sampled is below zero, zero is used (a pure ACO approach); if it is above one, one 
	 * is used (a pure Metropolis approach). If several solutions are created from the same values 
	 * generated using the ACO approach (see {@link #groupNumber} and {@link #groupPercent}), a 
	 * single ratio is used for the entire group.
	 */
	private ContProbDist independentRatio;
	
	/**
	 * The probability of each value of a generated solution to be re-sampled from a uniform 
	 * distribution spanning the entire valid range of the variable. Analogous to the mutation
	 * probability in Genetic Algorithms.
	 */
	private double uniformProb;
	
	/**
	 * The number of solutions to be generated from the same set of independently generated values
	 * using the ACO approach. 1 if each solution should be generated independently from the rest.
	 * The number of solutions in each group is still bounded by the number of solutions requested.
	 * Use a negative number if {@link #groupPercent} should be used instead. 
	 */
	private int groupNumber;
	
	/**
	 * The number of solutions to be generated from the same set of independently generated values
	 * (those using the ACO approach) as a ratio of the number of remaining decision variables to 
	 * be generated. The remaining variables are those to be generated using the Metropolis 
	 * approach. The number of solutions in each group is still bounded by the number of solutions 
	 * requested. Must be greater than zero.
	 */
	private double groupPercent;
	
	/**
	 * Determines the method used to compute the bandwidth matrix of the kernel density 
	 * distribution of the current population. The available methods are given by the following 
	 * class constants: <ul>
	 * <li>{@link #BANDWIDTH_SILVERMAN_UNIVAR}
	 * <li>{@link #BANDWIDTH_SCOTT}
	 * <li>{@link #BANDWIDTH_SILVERMAN} </ul>
	 */
	private int bandwidthMethod;
	
	/**
	 * Multiplier of the bandwidth values. After applying the method determined by 
	 * {@link #bandwidthMethod}, multiplies the resulting bandwidth values. Must be greater than
	 * zero.
	 */
	private double bandwidthMult;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new instance of the Metropolis-Ant Colony Optimization generator with the default 
	 * parameters
	 */
	public MetroACO()
	{
		discVars			= new ArrayList<DiscVar>();
		contVars			= new ArrayList<ContVar>();
		greed				= DEF_GREED;
		independentRatio	= DEF_INDEPENDENT_RATIO;
		uniformProb			= DEF_UNIFORM_PROB;
		groupNumber			= DEF_GROUP_NUMBER;
		groupPercent		= DEF_GROUP_PERCENT;
		bandwidthMethod		= DEF_BANDWIDTH_METHOD;
		bandwidthMult		= DEF_BANDWIDTH_MULT;
	}
	
	/**
	 * Creates a new instance of the Metropolis-Ant Colony Optimization generator
	 * @param greed 			{@link #greed}
	 * @param independentRatio	{@link #independentRatio}
	 * @param uniformProb		{@link #uniformProb}
	 * @param groupNumber		{@link #groupNumber}
	 * @param groupPercent		{@link #groupPercent}
	 * @param bandwidthMethod	{@link #bandwidthMethod}
	 * @param bandwidthMult		{@link #bandwidthMult}
	 */
	public MetroACO(double greed, ContProbDist independentRatio, double uniformProb, 
			int groupNumber, double groupPercent, int bandwidthMethod, double bandwidthMult)
	{
		discVars				= new ArrayList<DiscVar>();
		contVars				= new ArrayList<ContVar>();
		this.greed				= greed;
		this.independentRatio	= independentRatio;
		this.uniformProb		= uniformProb;
		this.groupNumber		= groupNumber;
		this.groupPercent		= groupPercent;
		this.bandwidthMethod	= bandwidthMethod;
		this.bandwidthMult		= bandwidthMult;
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
		this.greed = greed < -1.0 ? -1.0 : greed;
		this.greed = greed >  1.0 ?  1.0 : this.greed;
	}

	/**
	 * @return {@link #independentRatio}
	 */
	public ContProbDist getIndependentRatio() 
	{
		return independentRatio;
	}

	/**
	 * @param independentRatio {@link #independentRatio}
	 */
	public void setIndependentRatio(ContProbDist independentRatio) 
	{
		this.independentRatio = independentRatio;
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

	/**
	 * @return {@link #groupNumber}
	 */
	public int getGroupNumber() 
	{
		return groupNumber;
	}

	/**
	 * @param groupNumber {@link #groupNumber}
	 */
	public void setGroupNumber(int groupNumber) 
	{
		this.groupNumber = groupNumber;
	}

	/**
	 * @return {@link #groupPercent}
	 */
	public double getGroupPercent() 
	{
		return groupPercent;
	}

	/**
	 * @param groupPercent {@link #groupPercent}
	 */
	public void setGroupPercent(double groupPercent) 
	{
		this.groupPercent = groupPercent >= 0.0 ? 0.0 : groupPercent;
	}

	/**
	 * @return {@link #bandwidthMethod}
	 */
	public int getBandwidthMethod() 
	{
		return bandwidthMethod;
	}

	/**
	 * @param bandwidthMethod {@link #bandwidthMethod}
	 */
	public void setBandwidthMethod(int bandwidthMethod) 
	{
		this.bandwidthMethod = bandwidthMethod;
	}

	/**
	 * @return {@link #bandwidthMult}
	 */
	public double getBandwidthMult() 
	{
		return bandwidthMult;
	}

	/**
	 * @param bandwidthMult {@link #bandwidthMult}
	 */
	public void setBandwidthMult(double bandwidthMult) 
	{
		this.bandwidthMult = bandwidthMult <= 0.0 ? DEF_BANDWIDTH_MULT : bandwidthMult;
	}

	@Override
	public String getParamSummary() 
	{
		String bwMethod = bandwidthMethod == BANDWIDTH_SILVERMAN_UNIVAR ? 
				PARAM_B_W_SILVERMAN_UNIVAR : (bandwidthMethod == BANDWIDTH_SCOTT ? 
				PARAM_B_W_SCOTT : PARAM_B_W_SILVERMAN);
		
		String line	= "";
		line		+=			PARAM_GREED				+ " = " + greed;
		line		+= "; " +	PARAM_INDEPENDENT_RATIO + " = " + independentRatio;
		line		+= "; " +	PARAM_UNIFORM_PROB		+ " = " + uniformProb;
		line		+= "; " +	PARAM_GROUP_NUMBER		+ " = " + groupNumber;
		line		+= "; " +	PARAM_GROUP_PERCENT		+ " = " + groupPercent;
		line		+= "; " +	PARAM_BANDWIDTH_METHOD	+ " = " + bwMethod;
		line		+= "; " +	PARAM_BANDWIDTH_MULT	+ " = " + bandwidthMult;
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
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
	}
	
	@Override
	public ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population, int number)
	{
		int discVarCount		= discVars.size();
		int contVarCount		= contVars.size();
		int varCount			= discVarCount + contVarCount;
		
		// Compute weights
		double[] weights		= computeWeights(population);
		
		// Create probability distributions
		ArrayList<DiscProbDist> discDists	= null;
		MultiVarKernelDensity contDist		= null;
		if (discVarCount > 0) // Could not be necessary
			discDists						= createDiscreteDistributions(discVarCount, 
																			population, weights);
		if (contVarCount > 0)
			contDist						= createContinuousDistribution(contVarCount, 
																			population, weights);
		
		// Generate groups of solutions
		ArrayList<SolutionRoot> roots		= new ArrayList<>();
		while (roots.size() < number)
		{
			// Determine ACO/Metropolis ratio
			double independent	= independentRatio.sample();
			independent			= independent < 0.0 ? 0.0 : independent;
			independent			= independent > 1.0 ? 1.0 : independent;
			int acoCount		= (int)(independent*varCount);
			int metroCount		= varCount - acoCount;
			
			// Determine group size
			int groupSize		= groupNumber > 0 ? groupNumber : (int)(metroCount*groupPercent);
			int maxGen			= number - roots.size();
			groupSize			= groupSize < 1 		? 1 		: groupSize;
			groupSize			= groupSize > maxGen 	? maxGen	: groupSize;
			
			// Initialize the values for the base solution
			ArrayList<Integer> discVals	= null;
			ArrayList<Double> contVals	= null;
			if (discVarCount > 0)
				discVals				= new ArrayList<>();
			if (contVarCount > 0)
				contVals				= new ArrayList<>();
			ArrayList<Integer> acoDims	= new ArrayList<>();
			
			// ACO values generation
			if (acoCount > 0)
				generateACOValues(discVarCount, varCount, discDists, contDist, acoCount, discVals, 
									contVals, acoDims);
			
			// Metropolis values generation
			if (metroCount > 0)
			{
				if (acoCount == 0)
					roots.addAll(metropolisGenerate(discVarCount, contVarCount, varCount, 
													population, weights, contDist, groupSize));
				else
					roots.addAll(metroACOGenerate(discVarCount, contVarCount, varCount, population, 
							weights, contDist, acoCount, groupSize, discVals, contVals, acoDims));
			}
			else
				roots.add(new SolutionRoot(discVals, null));
		}
		
		// Apply uniform re-sampling
		if (uniformProb > 0.0)
			for (SolutionRoot root : roots)
				for (int i = 0; i < varCount; i++)
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
							int index	= i - discVarCount;
							ContVar var	= contVars.get(index);
							double val	= Uniform.sample(var.getMin(), var.getMax());
							root.getContValues().set(index, val);
						}
					}
		
		return roots;
	}

	/**
	 * Assigns a weight for each solution in the provided list according to the rank and the value
	 * of {@link #greed}
	 * @param solutions The list of solutions to assign the weights to
	 */
	private double[] computeWeights(ArrayList<Solution> solutions)
	{		
		// Obtain q
		greed				= greed >  1.0 ?  1.0 : greed;
		greed				= greed < -1.0 ? -1.0 : greed;		
		double temp			= Math.pow(1 - Math.abs(greed), GREED_TO_Q_POW);
		double q			= MIN_Q + (MAX_Q - MIN_Q)*temp;
		
		// Assign weights
		int solCount		= solutions.size();
		double[] weights	= new double[solCount];
		double stDev		= q*solCount;
		for (int i = 0; i < solCount; i++)
			weights[i]		= Normal.computepdf(1.0, stDev, i + 1);
		return weights;
	}

	/**
	 * Creates a series of multinomial distributions for each of the discrete variables of the 
	 * provided solutions (taking their weights into account)
	 * @param discVarCount	The number of discrete decision variables in the solutions
	 * @param solutions		The list of solutions to populate the distributions
	 */
	private ArrayList<DiscProbDist> createDiscreteDistributions(int discVarCount,
												ArrayList<Solution> solutions, double[] weights) 
	{
		ArrayList<ContSeries> values 	= new ArrayList<>();
		for (int i = 0; i < discVarCount; i++)
			values.add(new ContSeries(true));
		for (int i = 0; i < solutions.size(); i++)
		{
			Solution solution			= solutions.get(i);
			ArrayList<Integer> solVal	= solution.getDiscValues();
			for (int j = 0; j < discVarCount; j++)
				values.get(j).addValue(solVal.get(j), weights[i]);
		}
		ArrayList<DiscProbDist> discDists = new ArrayList<>();
		for (int j = 0; j < discVarCount; j++)
			discDists.add(new DiscProbDist(values.get(j)));
		return discDists;
	}
	
	/**
	 * Populates the continuous kernel density distribution with the values in the solutions 
	 * (taking their weights into account) and sets the bandwidth matrix
	 * @param contVarCount	The number of continuous decision variables in the solutions
	 * @param solutions		The list of solutions to populate the distribution
	 * @param contDist		The distribution to populate
	 */
	private MultiVarKernelDensity createContinuousDistribution(int contVarCount,
							ArrayList<Solution> solutions, double[] weights)
	{
		MultiVarKernelDensity contDist	= new MultiVarKernelDensity();
		contDist.setWeighted(true);
		for (int i = 0; i < solutions.size(); i++)
		{
			Solution solution			= solutions.get(i);
			ArrayList<Double> values	= solution.getContValues();
			Sample sample				= new Sample(weights[i], values);
			contDist.addSample(sample);
		}
		switch (bandwidthMethod)
		{
			case BANDWIDTH_SILVERMAN_UNIVAR:	contDist.silvermanUnivarDiagBW();	break;
			case BANDWIDTH_SCOTT:				contDist.scottDiagBW(); 				break;
			case BANDWIDTH_SILVERMAN:			contDist.silvermanDiagBW();			break;
			default:							contDist.silvermanDiagBW();			break;
		}
		if (bandwidthMult != 1.0)
			contDist.setBandwidth(MatUtil.multiply(contDist.getBandwidth(), bandwidthMult));
		return contDist;
	}
	
	/**
	 * Generates the values of a new solution using the independent marginalized distributions.
	 * @param discVarCount	The number of discrete variables
	 * @param varCount		The total number of variables
	 * @param discDists		{@link #createDiscreteDistributions}
	 * @param contDist		{@link #createContinuousDistribution}
	 * @param acoCount		The number of variables to sample using the ACO approach
	 * @param discVals		The vector to write the generated values to for the discrete variables. 
	 * Negative values are written if the values are to be generated using the Metropolis approach.
	 * @param contVals		The vector to write the generated values to for the continuous
	 * variables. {@link java.lang.Double#NaN} is written if the values are to be generated using 
	 * the Metropolis approach.
	 * @param acoDims		A list of the indices of the variables that are to be generated using 
	 * the ACO approach. The indices assume the continuous variable list is concatenated to the
	 * discrete variable list.
	 */
	private void generateACOValues(int discVarCount, int varCount, 
			ArrayList<DiscProbDist> discDists, MultiVarKernelDensity contDist,	int acoCount, 
			ArrayList<Integer> discVals, ArrayList<Double> contVals, ArrayList<Integer> acoDims)
	{
		// Select variables to generate using ACO
		for (int i = 0; i < varCount; i++)
			acoDims.add(i);
		Collections.shuffle(acoDims);
		for (int i = acoCount; i < varCount; i++)
			acoDims.remove(acoCount);
		Collections.sort(acoDims);
		
		// Sample values from independent distributions
		int dimIndex	= 0;
		for (int i = 0; i < varCount; i++)
		{
			if (i < discVarCount)
			{
				if (dimIndex < acoDims.size() && i == acoDims.get(dimIndex))
				{
					int value	= discDists.get(i).sample();
					discVals.add(discVars.get(i).validate(value));
					dimIndex++;
				}
				else
					discVals.add(-1);
			}
			else
			{
				if (dimIndex < acoDims.size() && i == acoDims.get(dimIndex))
				{
					KernelDensity marginal	= contDist.marginalize(i - discVarCount);
					double value			= marginal.sample();
					contVals.add(contVars.get(i - discVarCount).validate(value));
					dimIndex++;
				}
				else
					contVals.add(Double.NaN);
			}
		}
	}

	/**
	 * Generates a group of solution roots using the Metropolis approach
	 * @param discVarCount	The number of discrete variables
	 * @param contVarCount	The number of continuous variables
	 * @param varCount		The total number of variables
	 * @param solList		The list of solutions in the current population
	 * @param weights		The list of weight values for each solution
	 * @param contDist		{@link #createContinuousDistribution} 
	 * @param groupSize		The number of solution roots to generate
	 * {@link java.lang.Double#NaN} is written if the values are to be generated using the 
	 * Metropolis approach.
	 */
	private ArrayList<SolutionRoot> metropolisGenerate(int discVarCount, int contVarCount,
			int varCount, ArrayList<Solution> solList, double[] weights, 
			MultiVarKernelDensity contDist, int groupSize) 
	{
		ArrayList<SolutionRoot> roots	= new ArrayList<>();
		
		// Select root samples
		ContSeries selector				= new ContSeries(true);
		for (int j = 0; j < solList.size(); j++)
			selector.addValue(weights[j]);
		ArrayList<Double> solutionIDs	= selector.sampleMultiple(groupSize, true, true);
		
		// Generate values
		for (int j = 0; j < groupSize; j++)
		{
			Solution solution				= solList.get((int)(double)solutionIDs.get(j));
			ArrayList<Integer> discValsJ	= discVarCount > 0 ? solution.getDiscValues() : null;
			ArrayList<Double> contValsJ		= contVarCount > 0 ? new ArrayList<>() : null;
			if (contVarCount > 0)
			{
				ArrayList<Double> contSol	= solution.getContValues();
				double[] contSolA			= new double[contVarCount];
				for (int i = 0; i < contVarCount; i++)
					contSolA[i]				= contSol.get(i);
				double[][] bwMat			= contDist.getBandwidth();
				MultiVarNormal kernel		= new MultiVarNormal(contSolA, bwMat);
				double[] sample				= kernel.sample();
				for (int i = discVarCount; i < varCount; i++)
				{
					int index			= i - discVarCount;
					contValsJ.add(contVars.get(index).validate(sample[index]));
				}
			}
			roots.add(new SolutionRoot(discValsJ, contValsJ));
		}
		return roots;
	}

	/**
	 * Generates a group of solution roots using the Metropolis approach based on the ACO values
	 * provided
	 * @param discVarCount	The number of discrete variables
	 * @param contVarCount	The number of continuous variables
	 * @param varCount		The total number of variables
	 * @param solList		The list of solutions in the current population
	 * @param weights		The list of weight values for each solution
	 * @param contDist		{@link #createContinuousDistribution} 
	 * @param acoCount		The number of variables that were sampled using the ACO approach
	 * @param groupSize		The number of solution roots to generate
	 * @param discVals		The vector with the ACO-generated discrete values. Negative values are 
	 * written if the values are to be generated using the Metropolis approach.
	 * @param contVals		The vector with the ACO-generated continuous values. 
	 * {@link java.lang.Double#NaN} is written if the values are to be generated using the 
	 * Metropolis approach.
	 * @param acoDims		A list of the indices of the variables that were generated using the 
	 * ACO approach. The indices assume the continuous variable list is concatenated to the 
	 * discrete variable list.
	 */
	private ArrayList<SolutionRoot> metroACOGenerate(int discVarCount, int contVarCount,
			int varCount, ArrayList<Solution> solList, double[] weights, 
			MultiVarKernelDensity contDist, int acoCount, int groupSize, 
			ArrayList<Integer> discVals, ArrayList<Double> contVals, ArrayList<Integer> acoDims)
	{		
		// Prepare evaluation kernel
		MultiVarNormal kernel			= null;
		if (contVarCount > 0)
		{
			ArrayList<Integer> acoDimsT	= new ArrayList<>();
			double[] reducedMeans		= new double[acoCount];
			for (int i = 0; i < acoCount; i++)
			{
				int dimIndex			= acoDims.get(i);
				acoDimsT.add(dimIndex);
				reducedMeans[i]			= contVals.get(dimIndex);
			}
			double[][] reducedBWMat		= contDist.marginalize(acoDimsT).getBandwidth();
			kernel						= new MultiVarNormal(reducedMeans, reducedBWMat);
		}
		
		// Compute conditional weights
		ContSeries selector				= new ContSeries(true); 
		for (int j = 0; j < solList.size(); j++)
		{
			Solution solution			= solList.get(j);
			double weight				= 1.0;
			if (discVarCount > 0)
			{
				double resemblance		= 0;
				for (int i = 0; i < discVarCount; i++)
				{
					int value			= discVals.get(i);
					if (value >= 0)
					{
						int solVal		= solution.getDiscValues().get(i);
						DiscVar variable = discVars.get(i);
						if (variable.isScalar())
							resemblance	+= 1.0 - ((double)Math.abs(value - solVal))
														/(variable.getRange() - 1);
						else
							resemblance	+= value == solVal ? 1.0 : 0.0;
					}
				}
				weight					= resemblance/acoCount;
			}
			if (contVarCount > 0)
			{
				double[] reduced		= new double[acoCount];
				for (int i = 0; i < acoCount; i++)
					reduced[i]			= solution.getContValues().get(acoDims.get(i));
				double pdf				= kernel.getpdf(reduced);
				weight					*= Double.isNaN(pdf) ? 0.0 : pdf;
			}
			weight						*= weights[j];
			selector.addValue(j, weight);
		}
		
		// Generate values and create solution roots
		ArrayList<SolutionRoot> roots	= new ArrayList<>();
		ArrayList<Double> solutionIDs	= selector.sampleMultiple(groupSize, true, true);
		for (int j = 0; j < groupSize; j++)
		{
			Solution solution			= solList.get((int)(double)solutionIDs.get(j));
			ArrayList<Integer> discSol	= solution.getDiscValues();
			double[] sample				= new double[contVarCount];
			if (contVarCount > 0)
			{
				ArrayList<Double> contSol	= solution.getContValues();
				double[] contSolA			= new double[contVarCount];
				for (int i = 0; i < contVarCount; i++)
					contSolA[i]				= contSol.get(i);
				double[][] bwMat			= contDist.getBandwidth();
				kernel						= new MultiVarNormal(contSolA, bwMat);
				sample						= kernel.sample();
			}
			ArrayList<Integer> discValsJ	= discVarCount > 0 ? new ArrayList<>() : null;
			ArrayList<Double> contValsJ		= contVarCount > 0 ? new ArrayList<>() : null;
			for (int i = 0; i < varCount; i++)
			{
				if (i < discVarCount)
				{
					int base			= discVals.get(i);
					discValsJ.add(base < 0 ? discSol.get(i) : base);
				}
				else
				{
					int index			= i - discVarCount;
					double base			= contVals.get(index);
					contValsJ.add(Double.isNaN(base) ? contVars.get(index).validate(sample[index])
														: base);
				}
			}
			roots.add(new SolutionRoot(discValsJ, contValsJ));
		}
		return roots;
	}

}