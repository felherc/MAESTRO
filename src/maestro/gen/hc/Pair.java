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

import maestro.ContVar;
import maestro.DiscVar;
import maestro.gen.hc.Pair;
import maestro.gen.hc.SolutionRank;
import maestro.solution.SolutionRoot;
import probDist.DiscProbDist;
import probDist.Normal;

/**
 * This class represents a pair of solutions that define a gradient within a search space in an
 * optimization problem. The gradient is used to generate new candidate solutions.
 * @author Felipe Hernández
 */
public class Pair implements Comparable<Pair>
{
	
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
	 * The high-rank solution in the pair
	 */
	private SolutionRank highSolution;
	
	/**
	 * The low-rank solution in the pair
	 */
	private SolutionRank lowSolution;
	
	/**
	 * The ratio between the total difference of the solutions and the difference in rank. The 
	 * total difference is the sum of the normalized differences for each variable. Pairs with 
	 * larger gradients are preferred to generate new solutions.
	 */
	private double gradient;
	
	/**
	 * Integer accumulation slot associated with the pair
	 */
	private int hits;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param discVars The list of the discrete variables of the problem
	 * @param contVars The list of the continuous variables of the problem
	 * @param solution1 The first solution of the pair
	 * @param solution2 The second solution of the pair
	 */
	public Pair(ArrayList<DiscVar> discVars, ArrayList<ContVar> contVars, 
				SolutionRank solution1, SolutionRank solution2)
	{
		this.discVars = discVars;
		this.contVars = contVars;
		if(solution1.rank < solution2.rank)
		{
			highSolution = solution1;
			lowSolution = solution2;
		}
		else
		{
			highSolution = solution2;
			lowSolution = solution1;
		}
		computeGradient();
		hits = 0;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * Computes the gradient of the pair
	 */
	private void computeGradient()
	{
		double totalDifference = 0;
		for(int i = 0 ; i < discVars.size() ; i++)
		{
			DiscVar var = discVars.get(i);
			if(var.isScalar())
			{
				double valueH = highSolution.solution.getDiscValues().get(i);
				double valueL = lowSolution.solution.getDiscValues().get(i);
				totalDifference += Math.abs((valueH - valueL)/(double)(var.getCount() - 1));
			}
			else
				totalDifference += highSolution.solution.getDiscValues().get(i) == 
									lowSolution.solution.getDiscValues().get(i) ? 0.0 : 1.0;
		}
		for(int i = 0 ; i < contVars.size() ; i++)
		{
			ContVar var = contVars.get(i);
			double valueH = highSolution.solution.getContValues().get(i);
			double valueL = lowSolution.solution.getContValues().get(i);
			totalDifference += Math.abs((valueH - valueL)/(var.getMax() - var.getMin()));
		}
		gradient = totalDifference/(double)(lowSolution.rank - highSolution.rank);
	}
	
	/**
	 * @return The high-rank solution in the pair
	 */
	public SolutionRank getHighSolution() 
	{
		return highSolution;
	}

	/**
	 * @return The low-rank solution in the pair
	 */
	public SolutionRank getLowSolution() 
	{
		return lowSolution;
	}

	/**
	 * @return The ratio between the total difference of the solutions and the difference in rank. 
	 * The total difference is the sum of the normalized differences for each variable. Pairs with 
	 * larger gradients are preferred to generate new solutions.
	 */
	public double getGradient() 
	{
		return gradient;
	}
	
	/**
	 * @return Integer accumulation slot associated with the pair
	 */
	public int getHits() 
	{
		return hits;
	}

	/**
	 * @param hits The number of hits to accumulate
	 */
	public void addHits(int hits) 
	{
		this.hits += hits;
	}
	
	/**
	 * Resets the integer accumulation slot associated with the pair
	 */
	public void resetHits()
	{
		hits = 0;
	}
	
	/**
	 * Generates new solutions roots from the gradient created by the two solutions in the pair
	 * @param count The number of solution roots to generate
	 * @param extent The extent parameter of the search range for new solutions. New solutions are 
	 * generated using a normal distribution to randomly sample each of its new values. The mean of 
	 * the distribution is created by projecting the gradient of the high rank parent (<i>p1</i>) 
	 * and the low rank parent (<i>p2</i>) for each variable <i>i</i> using the extent as a 
	 * percentage of the difference: <i>mean = p1i + extent*(p1i - p2i)</i>.
	 * @param amplitude The amplitude parameter of the search range for new solutions. New 
	 * solutions are generated using a normal distribution to randomly sample each of its new 
	 * values. The amplitude is the standard deviation of the distribution as a percentage of 
	 * <i>p1i - p2i</i>.
	 * @return The solutions roots generated
	 */
	public ArrayList<SolutionRoot> generateSolutions(int count, double extent, double amplitude)
	{
		// Create distributions
		ArrayList<DiscProbDist> discDistributions = new ArrayList<DiscProbDist>();
		ArrayList<Normal> contDistributions = new ArrayList<Normal>();
		for(int i = 0 ; i < discVars.size() ; i++)
		{
			DiscVar var = discVars.get(i);
			if(var.isScalar())
			{
				double high = highSolution.solution.getDiscValues().get(i);
				double low = lowSolution.solution.getDiscValues().get(i);
				double mean = high + extent*(high - low);
				double stDev = amplitude*(high - low);
				Normal contDist = new Normal(mean, stDev);
				discDistributions.add(new DiscProbDist(var.getMin(), var.getMax(), contDist));
			}
			else
			{
				int high = highSolution.solution.getDiscValues().get(i);
				int low = lowSolution.solution.getDiscValues().get(i);
				if(high == low)
					discDistributions.add(new DiscProbDist(high, high));
				else
					discDistributions.add(new DiscProbDist(var.getMin(), var.getMax()));
			}
		}
		for(int i = 0 ; i < contVars.size() ; i++)
		{
			double high = highSolution.solution.getContValues().get(i);
			double low = lowSolution.solution.getContValues().get(i);
			double mean = high + extent*(high - low);
			double stDev = amplitude*(high - low);
			contDistributions.add(new Normal(mean, stDev));
		}
		
		// Generate solutions roots
		ArrayList<SolutionRoot> roots = new ArrayList<SolutionRoot>();
		for(int j = 0 ; j < count ; j++)
		{
			ArrayList<Integer> discValues = new ArrayList<Integer>();
			ArrayList<Double> contValues = new ArrayList<Double>();
			for(int i = 0 ; i < discVars.size() ; i++)
				discValues.add(discDistributions.get(i).sample());
			for(int i = 0 ; i < contVars.size() ; i++)
			{
				ContVar var = contVars.get(i);
				double value = contDistributions.get(i).sample();
				value = value < var.getMin() ? var.getMin() : value;
				value = value > var.getMax() ? var.getMax() : value;
				contValues.add(value);
			}
			roots.add(new SolutionRoot(discValues, contValues));
		}

		return roots;
	}
	
	@Override
	public int compareTo(Pair other)
	{
		Double thisGrad = this.gradient;
		Double otherGrad = other.getGradient();
		return thisGrad.compareTo(otherGrad);
	}
	
}
