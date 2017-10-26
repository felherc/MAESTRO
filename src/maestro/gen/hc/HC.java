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

package maestro.gen.hc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.TreeSet;

import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.Generator;
import maestro.gen.hc.Pair;
import maestro.gen.hc.SolutionRank;
import maestro.solution.Solution;
import maestro.solution.SolutionRoot;
import utilities.Utilities;

/**
 * This generator allows to create new solutions to a problem with discrete and continuous 
 * variables using a sorted base population of solutions and a Hill-Climbing algorithm
 * @author Felipe Hernández
 */
public class HC implements Generator
{
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator method
	 */
	public final static String ID = "Hill-Climbing";
	
	/**
	 * The short version of the identifier of the generator method
	 */
	public final static String SHORT_ID = "HC";
	
	/**
	 * The <code>qSolutions</code> attribute printing name
	 */
	public final static String PARAM_Q_SOLUTIONS = "q solutions = ";
	
	/**
	 * The <code>qPairs</code> attribute printing name
	 */
	public final static String PARAM_Q_PAIRS = "q pairs = ";
	
	/**
	 * The <code>truncSolutions</code> attribute printing name
	 */
	public final static String PARAM_TRUNC_SOLUTIONS = "trunc solutions = ";
	
	/**
	 * The <code>truncPairs</code> attribute printing name
	 */
	public final static String PARAM_TRUNC_PAIRS = "trunc pairs = ";
	
	/**
	 * The <code>selectionMethodSolutions</code> attribute printing name
	 */
	public final static String PARAM_SELECTION_METHOD_SOLUTIONS = "Solution selection method = ";
	
	/**
	 * The <code>selectionMethodPairs</code> attribute printing name
	 */
	public final static String PARAM_SELECTION_METHOD_PAIRS = "Pair selection method = ";
	
	/**
	 * The <code>extent</code> attribute printing name
	 */
	public final static String PARAM_EXTENT = "extent = ";
	
	/**
	 * The <code>amplitude</code> attribute printing name
	 */
	public final static String PARAM_AMPLITUDE = "amplitude = ";
	
	/**
	 * Fitness proportionate selection or Roulette-wheel selection printing name
	 */
	public final static String SEL_ROULETTE_PRINT = "roulette-wheel";
	
	/**
	 * Stochastic universal sampling selection printing name
	 */
	public final static String SEL_SUS_PRINT = "stochastic universal sampling";
	
	/**
	 * Tournament selection printing name
	 */
	public final static String SEL_TOURNAMENT_PRINT = "tournament";
	
	/**
	 * Index of the Fitness proportionate selection or Roulette-wheel selection strategy
	 */
	public final static int SELECTION_ROULETTE = 0;
	
	/**
	 * Index of the Stochastic universal sampling selection strategy
	 */
	public final static int SELECTION_SUS = 1;
	
	/**
	 * Index of the Tournament selection strategy
	 */
	public final static int SELECTION_TOURNAMENT = 2;
	
	/**
	 * Default value for the <code>qSolutions</code> attribute
	 */
	public final static double D_Q_SOLUTIONS = 0.4; 
	// Moderate importance
	
	/**
	 * Default value for the <code>qPairs</code> attribute
	 */
	public final static double D_Q_PAIRS = 0.2; 
	// Low importance
	
	/**
	 * Default value for the <code>truncSolutions</code> attribute
	 */
	public final static double D_TRUNC_SOLUTIONS = 0.3; 
	// Insignificant - near sweet spot
	
	/**
	 * Default value for the <code>truncPairs</code> attribute
	 */
	public final static double D_TRUNC_PAIRS = 0.6; 
	// Insignificant - above 0.2 is better
	
	/**
	 * Default value for the <code>selectionMethodSolutions</code> attribute
	 */
	public final static int D_SELECTION_METHOD_SOLUTIONS = 0; 
	// High importance - Best: 0
	
	/**
	 * Default value for the <code>selectionMethodPairs</code> attribute
	 */
	public final static int D_SELECTION_METHOD_PAIRS = 0; 
	// Low importance - Best: 0
	
	/**
	 * Default value for the <code>extent</code> attribute
	 */
	public final static double D_EXTENT = 1.0; 
	// Highest importance - near sweet spot (10); smaller is better for more dimension
	
	/**
	 * Default value for the <code>amplitude</code> attribute
	 */
	public final static double D_AMPLITUDE = 1.25; 
	// Moderate importance - interaction with extent
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The list of the discrete variables of the problem
	 */
	private ArrayList<DiscVar> discVars;
	
	/**
	 * The list of the continuous variables of the problem
	 */
	private ArrayList<ContVar> contVars;
	
	/**
	 * This parameter defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>qSolutions</i> has to be 
	 * positive. A value of <i>0</i> means only the most fitted solution may be selected. A large 
	 * value of <i>qSolutions</i> means that every solution has the same probability of being 
	 * selected with little regard for its ranking. Use <code>Double.NaN</code> if an exact uniform 
	 * distribution should be used.
	 */
	private double qSolutions;
	
	/**
	 * This parameter defines how the solution pairs are assigned weights. The weights are used to 
	 * compute the likelihood of being selected. <i>qPairs</i> has to be positive. A value of 
	 * <i>0</i> means only the most fitted pair may be selected. A large value of <i>qPairs</i> 
	 * means that every pair has the same probability of being selected with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 */
	private double qPairs;
	
	/**
	 * The percentage of the number of solutions that are to be selected for constructing pairs. 
	 * <i>(1 - <code>truncSolutions</code>)*n</i> solutions are selected.
	 */
	private double truncSolutions;
	
	/**
	 * The percentage of the number of pairs that are to be selected for constructing new 
	 * solutions. The <i>(1 - <code>truncPairs</code>)*n</i> best pairs are selected.
	 */
	private double truncPairs;
	
	/**
	 * The index of the selection method for the solutions as defined by the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	private int selectionMethodSolutions;
	
	/**
	 * The index of the selection method for the pairs as defined by the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	private int selectionMethodPairs;
	
	/**
	 * The extent parameter of the search range for new solutions. New solutions are generated
	 * using a normal distribution to randomly sample each of its new values. The mean of the 
	 * distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) and 
	 * the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a percentage
	 * of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 */
	private double extent;
	
	/**
	 * The amplitude parameter of the search range for new solutions. New solutions are generated
	 * using a normal distribution to randomly sample each of its new values. The amplitude is the 
	 * standard deviation of the distribution as a percentage of <i>p1i - p2i</i>.
	 */
	private double amplitude;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new instance of the Hill-Climbing generator with the default parameters
	 */
	public HC()
	{
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
		qSolutions = D_Q_SOLUTIONS;
		qPairs = D_Q_PAIRS;
		truncSolutions = D_TRUNC_SOLUTIONS;
		truncPairs = D_TRUNC_PAIRS;
		selectionMethodSolutions = D_SELECTION_METHOD_SOLUTIONS;
		selectionMethodPairs = D_SELECTION_METHOD_PAIRS;
		extent = D_EXTENT;
		amplitude = D_AMPLITUDE;
	}
	
	/**
	 * Creates a new instance of theHill-Climbing generator
	 * @param qSolutions This parameter defines how the solutions in the population are assigned 
	 * weights. The weights are used to compute the likelihood of being selected. <i>qSolutions</i> 
	 * has to be positive. A value of <i>0</i> means only the most fitted solution may be selected. 
	 * A large value of <i>qSolutions</i> means that every solution has the same probability of 
	 * being selected with little regard for its ranking. Use <code>Double.NaN</code> if an exact 
	 * uniform distribution should be used.
	 * @param qPairs This parameter defines how the solution pairs are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>qPairs</i> has to be 
	 * positive. A value of <i>0</i> means only the most fitted pair may be selected. A large value 
	 * of <i>qPairs</i> means that every pair has the same probability of being selected with 
	 * little regard for its ranking. Use <code>Double.NaN</code> if an exact uniform distribution 
	 * should be used.
	 * @param truncSolutions The percentage of the number of solutions that are to be selected for 
	 * constructing pairs. <i>(1 - <code>truncSolutions</code>)*n</i> solutions are selected.
	 * @param truncPairs The percentage of the number of pairs that are to be selected for 
	 * constructing new solutions. The <i>(1 - <code>truncPairs</code>)*n</i> best pairs are 
	 * selected.
	 * @param selectionMethodSolutions The index of the selection method for the solutions as 
	 * defined by the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 * @param slectionMethodPairs The index of the selection method for the pairs as defined by the 
	 * class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 * @param extent The extent parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The mean of 
	 * the distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) and 
	 * the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a percentage
	 * of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 * @param amplitude The amplitude parameter of the search range for new solutions. New 
	 * solutions are generated using a normal distribution to randomly sample each of its new 
	 * values. The amplitude is the standard deviation of the distribution as a percentage of 
	 * <i>p1i - p2i</i>.
	 */
	public HC(double qSolutions, double qPairs, double truncSolutions, double truncPairs, 
				int selectionMethodSolutions, int selectionMethodPairs, double extent, 
				double amplitude)
	{
		discVars = new ArrayList<DiscVar>();
		contVars = new ArrayList<ContVar>();
		if(Double.isNaN(qSolutions))
			this.qSolutions = qSolutions;
		else
			this.qSolutions = qSolutions < 0 ? 0 : qSolutions;
		if(Double.isNaN(qPairs))
			this.qPairs = qPairs;
		else
			this.qPairs = qPairs < 0 ? 0 : qPairs;
		this.truncSolutions = truncSolutions < 0 ? 0.0 : truncSolutions;
		this.truncSolutions = truncSolutions > 1 ? 1.0 : this.truncSolutions;
		this.truncPairs = truncPairs < 0 ? 0.0 : truncPairs;
		this.truncPairs = truncPairs > 1 ? 1.0 : this.truncPairs;
		this.selectionMethodSolutions = selectionMethodSolutions;
		this.selectionMethodPairs = selectionMethodPairs;
		this.extent = extent;
		this.amplitude = amplitude < 0 ? 0 : amplitude;
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
	 * @return This parameter defines how the solutions in the population are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>qSolutions</i> has to be 
	 * positive. A value of <i>0</i> means only the most fitted solution may be selected. A large 
	 * value of <i>qSolutions</i> means that every solution has the same probability of being 
	 * selected with little regard for its ranking. Use <code>Double.NaN</code> if an exact uniform 
	 * distribution should be used.
	 */
	public double getqSolutions() 
	{
		return qSolutions;
	}

	/**
	 * @param qSolutions This parameter defines how the solutions in the population are assigned 
	 * weights. The weights are used to compute the likelihood of being selected. <i>qSolutions</i> 
	 * has to be positive. A value of <i>0</i> means only the most fitted solution may be selected. 
	 * A large value of <i>qSolutions</i> means that every solution has the same probability of 
	 * being selected with little regard for its ranking. Use <code>Double.NaN</code> if an exact 
	 * uniform distribution should be used.
	 */
	public void setqSolutions(double qSolutions) 
	{
		if(Double.isNaN(qSolutions))
			this.qSolutions = qSolutions;
		else
			this.qSolutions = qSolutions < 0 ? 0 : qSolutions;
	}

	/**
	 * @return This parameter defines how the solution pairs are assigned weights. The weights are 
	 * used to compute the likelihood of being selected. <i>qPairs</i> has to be positive. A value 
	 * of <i>0</i> means only the most fitted pair may be selected. A large value of <i>qPairs</i> 
	 * means that every pair has the same probability of being selected with little regard for its 
	 * ranking. Use <code>Double.NaN</code> if an exact uniform distribution should be used.
	 */
	public double getqPairs() 
	{
		return qPairs;
	}

	/**
	 * @param qPairs This parameter defines how the solution pairs are assigned weights. The 
	 * weights are used to compute the likelihood of being selected. <i>qPairs</i> has to be 
	 * positive. A value of <i>0</i> means only the most fitted pair may be selected. A large value 
	 * of <i>qPairs</i> means that every pair has the same probability of being selected with 
	 * little regard for its ranking. Use <code>Double.NaN</code> if an exact uniform distribution 
	 * should be used.
	 */
	public void setqPairs(double qPairs) 
	{
		if(Double.isNaN(qPairs))
			this.qPairs = qPairs;
		else
			this.qPairs = qPairs < 0 ? 0 : qPairs;
	}

	/**
	 * @return The percentage of the number of solutions that are to be selected for constructing 
	 * pairs. <i>(1 - <code>truncSolutions</code>)*n</i> solutions are selected.
	 */
	public double getTruncSolutions() 
	{
		return truncSolutions;
	}

	/**
	 * @param truncSolutions The percentage of the number of solutions that are to be selected for 
	 * constructing pairs. <i>(1 - <code>truncSolutions</code>)*n</i> solutions are selected.
	 */
	public void setTruncSolutions(double truncSolutions) 
	{
		this.truncSolutions = truncSolutions < 0 ? 0.0 : truncSolutions;
		this.truncSolutions = truncSolutions > 1 ? 1.0 : this.truncSolutions;
	}

	/**
	 * @return The percentage of the number of pairs that are to be selected for constructing new 
	 * solutions. The <i>(1 - <code>truncPairs</code>)*n</i> best pairs are selected.
	 */
	public double getTruncPairs() 
	{
		return truncPairs;
	}

	/**
	 * @param truncPairs The percentage of the number of pairs that are to be selected for 
	 * constructing new solutions. The <i>(1 - <code>truncPairs</code>)*n</i> best pairs are 
	 * selected.
	 */
	public void setTruncPairs(double truncPairs) 
	{
		this.truncPairs = truncPairs < 0 ? 0.0 : truncPairs;
		this.truncPairs = truncPairs > 1 ? 1.0 : this.truncPairs;
	}
	
	/**
	 * @return The index of the selection method for the solutions as defined by the class 
	 * constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	public int getSelectionMethodSolutions() 
	{
		return selectionMethodSolutions;
	}

	/**
	 * @param selectionMethodSolutions The index of the selection method for the solutions as 
	 * defined by the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	public void setSelectionMethodSolutions(int selectionMethodSolutions) 
	{
		this.selectionMethodSolutions = selectionMethodSolutions;
	}

	/**
	 * @return The index of the selection method for the pairs as defined by the class constants: 
	 * <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	public int getSelectionMethodPairs() 
	{
		return selectionMethodPairs;
	}

	/**
	 * @param selectionMethodPairs The index of the selection method for the pairs as defined by 
	 * the class constants: <ul>
	 * <li> <code>SELECTION_ROULETTE</code>: Fitness proportionate selection or Roulette-wheel 
	 * selection
	 * <li> <code>SELECTION_SUS</code>: Stochastic universal sampling selection
	 * <li> <code>SELECTION_TOURNAMENT</code>: Tournament selection </ul>
	 */
	public void setSelectionMethodPairs(int selectionMethodPairs) 
	{
		this.selectionMethodPairs = selectionMethodPairs;
	}

	/**
	 * @return The extent parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The mean of 
	 * the distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) and 
	 * the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a percentage
	 * of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 */
	public double getExtent() 
	{
		return extent;
	}

	/**
	 * @param extent The extent parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The mean of 
	 * the distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) and 
	 * the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a percentage
	 * of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 */
	public void setExtent(double extent) 
	{
		this.extent = extent;
	}

	/**
	 * @return The amplitude parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The 
	 * amplitude is the standard deviation of the distribution as a percentage of <i>p1i - p2i</i>.
	 */
	public double getAmplitude() 
	{
		return amplitude;
	}

	/**
	 * @param amplitude The amplitude parameter of the search range for new solutions. New 
	 * solutions are generated using a normal distribution to randomly sample each of its new 
	 * values. The amplitude is the standard deviation of the distribution as a percentage of 
	 * <i>p1i - p2i</i>.
	 */
	public void setAmplitude(double amplitude) 
	{
		this.amplitude = amplitude < 0 ? 0 : amplitude;
	}

	@Override
	public String getParamSummary() 
	{
		String line = PARAM_Q_SOLUTIONS + qSolutions;
		line += "; " + PARAM_Q_PAIRS + qPairs;
		line += "; " + PARAM_TRUNC_SOLUTIONS + truncSolutions;
		line += "; " + PARAM_TRUNC_PAIRS + truncPairs;
		String solSelMethod = selectionMethodSolutions == SELECTION_ROULETTE ? SEL_ROULETTE_PRINT :
				(selectionMethodSolutions == SELECTION_SUS ? SEL_SUS_PRINT :
				(selectionMethodSolutions == SELECTION_TOURNAMENT ? SEL_TOURNAMENT_PRINT : "N/A"));
		line += "; " + PARAM_SELECTION_METHOD_SOLUTIONS + solSelMethod;
		String pairSelMethod = selectionMethodPairs == SELECTION_ROULETTE ? SEL_ROULETTE_PRINT :
					(selectionMethodPairs == SELECTION_SUS ? SEL_SUS_PRINT :
					(selectionMethodPairs == SELECTION_TOURNAMENT ? SEL_TOURNAMENT_PRINT : "N/A"));
		line += "; " + PARAM_SELECTION_METHOD_PAIRS + pairSelMethod;
		line += "; " + PARAM_EXTENT + extent;
		line += "; " + PARAM_AMPLITUDE + amplitude;
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
	public synchronized ArrayList<SolutionRoot> generateSolutions(ArrayList<Solution> population, 
																	int number)
	{
		// Select solutions
		ArrayList<Double> weights = computeWeights(population.size());
		int selectCount = Math.max(2, (int)(truncSolutions*population.size()));
		ArrayList<Integer> selectIds = null;
		switch(selectionMethodSolutions)
		{
			case SELECTION_ROULETTE:	selectIds = selectRoulette(weights, selectCount, 
																Double.isNaN(qSolutions)); break;
			case SELECTION_SUS:			selectIds = selectSUS(weights, selectCount); break;
			case SELECTION_TOURNAMENT:	selectIds = selectTournament(weights, selectCount); break;
			default:					selectIds = selectRoulette(weights, selectCount, 
																Double.isNaN(qSolutions)); break;
		}
		
		// Create pairs
		TreeSet<Pair> allPairs = new TreeSet<Pair>();
		int maxCount = 0;
		for(int i = 0 ; i < selectCount ; i++)
			for(int j = i + 1 ; j < selectCount ; j++)
			{
				int id1 = selectIds.get(i);
				int id2 = selectIds.get(j);
				if(id1 != id2)
				{
					SolutionRank sol1 = new SolutionRank(population.get(id1), id1);
					SolutionRank sol2 = new SolutionRank(population.get(id2), id2);
					Pair pair = new Pair(discVars, contVars, sol1, sol2);
					allPairs.add(pair);
				}
				maxCount++;
			}
		int pairSelectCount = Math.max(1, (int)(truncPairs*maxCount));
		ArrayList<Pair> pairs = new ArrayList<Pair>();
		Iterator<Pair> iter = allPairs.descendingIterator();
		if(allPairs.size() == 0)
		{
			SolutionRank sol1 = new SolutionRank(population.get(0), 0);
			SolutionRank sol2 = new SolutionRank(population.get(1), 1);
			pairs.add(new Pair(discVars, contVars, sol1, sol2));
		}
		else
			for(int i = 0 ; i < pairSelectCount ; i++)
				if(iter.hasNext())
					pairs.add(iter.next());
		
		// Select pairs
		weights = computeWeights(pairs.size());
		selectIds = null;
		switch(selectionMethodSolutions)
		{
			case SELECTION_ROULETTE:	selectIds = selectRoulette(weights, number, 
																Double.isNaN(qSolutions)); break;
			case SELECTION_SUS:			selectIds = selectSUS(weights, number); break;
			case SELECTION_TOURNAMENT:	selectIds = selectTournament(weights, number); break;
			default:					selectIds = selectRoulette(weights, number, 
																Double.isNaN(qSolutions)); break;
		}
		for(Integer pairId : selectIds)
			pairs.get(pairId).addHits(1);
		
		// Generate solution roots
		ArrayList<SolutionRoot> roots = new ArrayList<SolutionRoot>();
		for(Pair pair : pairs)
			roots.addAll(pair.generateSolutions(pair.getHits(), extent, amplitude));
		return roots;
	}
	
	/**
	 * Computes the weighting factors according to the value of <i>q</i> and the size of the 
	 * population. These weights are used for computing the probability of selecting a given 
	 * solution as a parent to generate new offspring.
	 * @param popSize The size of the population
	 */
	private ArrayList<Double> computeWeights(int popSize)
	{
		ArrayList<Double> weights = null;
		
		// Compute weights
		if(Double.isNaN(qSolutions))
		{
			weights = new ArrayList<Double>();
			for(int i = 0 ; i < popSize ; i++)
				weights.add(1.0);
		}
		else if(qSolutions > 0)
		{
			weights = new ArrayList<Double>();
			double sum = 0;
			for(int i = 0 ; i < popSize ; i++)
			{
				int rank = i + 1;
				double weight = 1/(qSolutions*popSize*Math.sqrt(2*Math.PI))
								* Math.exp(-(rank - 1)*(rank - 1)/
										(2*qSolutions*qSolutions*popSize*popSize));
				weights.add(weight);
				sum += weight;
			}
			
			// Normalize
			for(int i = 0 ; i < popSize ; i++)
				weights.set(i, weights.get(i)/sum);
		}
		return weights;
	}
	
	/**
	 * Returns the indices of the candidates selected using the Fitness proportionate selection 
	 * method
	 * @param weights The list with the probability of choosing each candidate (the fitness values)
	 * @param number The number of selected candidates to generate
	 * @param uniform True if all the candidates have the same weight (fitness value)
	 * @return The indices of the candidates selected using the Fitness proportionate selection 
	 * method
	 */
	private ArrayList<Integer> selectRoulette(ArrayList<Double> weights, int number, 
												boolean uniform)
	{
		ArrayList<Integer> parentIds = new ArrayList<Integer>();
		for(int gen = 0 ; gen < number ; gen++)
		{
			if(uniform)
				parentIds.add(Utilities.uniformRandomSelect(weights.size()));
			else
			{
				double rand = 1 - Math.random();
				double accum = 0;
				int index = 0;
				boolean ok = false;
				while(!ok)
				{
					accum += weights.get(index);
					if(rand <= accum || index == weights.size() - 1)
						ok = true;
					else
						index++;
				}
				parentIds.add(index);
			}
		}
		return parentIds;
	}
	
	/**
	 * Returns the indices of the candidates selected using the Stochastic universal sampling 
	 * selection method
	 * @param weights The list with the probability of choosing each candidate (the fitness values)
	 * @param number The number of selected candidates to generate
	 * @return The indices of the candidates selected using the Stochastic universal sampling 
	 * selection method
	 */
	private ArrayList<Integer> selectSUS(ArrayList<Double> weights, int number)
	{
		ArrayList<Integer> parentIds = new ArrayList<Integer>();
		double delta = Math.random();
		double pointer = 0;
		double accum = 0;
		int index = -1;
		for(int gen = 0 ; gen < number ; gen++)
		{
			pointer += delta;
			if(pointer > 1)
			{
				pointer = pointer - 1;
				index = -1;
				accum = 0;
			}
			while(pointer >= accum)
			{
				index++;
				accum += weights.get(index);
			}
			parentIds.add(index);
		}
		Collections.shuffle(parentIds);
		return parentIds;
	}
	
	/**
	 * Returns the indices of the candidates selected using the Tournament selection method
	 * @param weights The list with the probability of choosing each candidate (the fitness values)
	 * @param number The number of selected candidates to generate
	 * @return The indices of the candidates selected using the Tournament selection method
	 */
	private ArrayList<Integer> selectTournament(ArrayList<Double> weights, int number)
	{
		ArrayList<Integer> parentIds = new ArrayList<Integer>();
		for(int gen = 0 ; gen < number ; gen++)
		{
			ArrayList<Integer> partial = new ArrayList<Integer>();
			for(int i = 0 ; i < weights.size() ; i++)
			{
				if(Double.isNaN(qSolutions))
					partial.add(Utilities.uniformRandomSelect(weights.size()));
				else
				{
					double rand = 1 - Math.random();
					double accum = 0;
					int index = 0;
					boolean ok = false;
					while(!ok)
					{
						accum += weights.get(index);
						if(rand <= accum || index == weights.size() - 1)
							ok = true;
						else
							index++;
					}
					partial.add(index);
				}
			}
			if(partial.size() > 1)
				parentIds.add(Collections.min(partial));
			else
				parentIds.add(partial.get(0));
		}
		return parentIds;
	}

}
