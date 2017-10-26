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

package maestro;

import java.util.ArrayList;
import java.util.Collections;

import probDist.ContProbDist;
import probDist.Uniform;
import utilities.stat.ContSeries;

/**
 * Stores information for a continuous variable of the problem, including its name and the minimum
 * and maximum possible values
 * @author Felipe Hernández
 */
public class ContVar 
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the variable
	 */
	private String name;
	
	
	/**
	 * The minimum value the variable can take
	 */
	private double min;
	
	/**
	 * The maximum value the variable can take
	 */
	private double max;
	
	/**
	 * The probability distribution that fits the values of the solutions in the current population 
	 * set. The values are those corresponding to this variable. The distribution is used to 
	 * compute the rarity of new solutions.
	 */
	private ContProbDist distribution;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes a new continuous variable
	 * @param name Identifier of the variable
	 * @param min The minimum value the variable can take
	 * @param max The maximum value the variable can take
	 */
	public ContVar(String name, double min, double max)
	{
		this.name = name;
		this.min = min <= max ? min : max;
		this.max = min <= max ? max : min;
		distribution = null;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * @return Identifier of the variable
	 */
	public String getName() 
	{
		return name;
	}

	/**
	 * @return The minimum value the variable can take
	 */
	public double getMin() 
	{
		return min;
	}

	/**
	 * @return The maximum value the variable can take
	 */
	public double getMax() 
	{
		return max;
	}
	
	/**
	 * @return The value range of the variable (maximum value - minimum value)
	 */
	public double getRange()
	{
		return max - min;
	}

	/**
	 * @return The probability distribution that fits the values of the solutions in the current 
	 * population set. The values are those corresponding to this variable. The distribution is 
	 * used to compute the rarity of new solutions.
	 */
	public ContProbDist getDistribution() 
	{
		return distribution;
	}
	
	/**
	 * Generates random values to initialize the population
	 * @param count The number of values (solutions) to generate
	 * @param uniform True if the Latin hypercube sampling (LHS) method should be used to guarantee
	 * a uniform sampling of the variable's range. False to use random sampling.
	 * @return The random values generated
	 */
	public ArrayList<Double> generateRandomValues(int count, boolean uniform)
	{
		ArrayList<Double> values = new ArrayList<Double>();
		if(uniform)
		{
			for(int i = 0 ; i < count ; i ++)
			{
				double lowerBound = min + (max - min)*i/count;
				double upperBound = min + (max - min)*(i + 1)/count;
				values.add(Uniform.sample(lowerBound, upperBound));
			}
			Collections.shuffle(values);
		}
		else
			for(int i = 0 ; i < count ; i ++)
				values.add(Uniform.sample(min, max));
		return values;
	}
	
	/**
	 * Validates that a candidate value is between {@link #min} and {@link #max}
	 * @param candidate The value to validate
	 * @return The validated value
	 */
	public double validate(double candidate)
	{
		double validated	= candidate < min ? min : candidate;
		validated			= validated > max ? max : validated;
		return validated;
	}
	
	/**
	 * Adjust a probability distribution to the values observed in order to be able to compute the
	 * rarity of new solutions.
	 * @param values The values of the solutions in the current population set that correspond to 
	 * the variable
	 */
	public void adjustDistribution(ContSeries values)
	{
		distribution = ContProbDist.fromValues(values);
	}
	
	/**
	 * Computes the rarity of the provided value. The rarity measures how different or rare is the 
	 * value in comparison with the values in the first partition of the current population. The 
	 * rarity is a value form 0 (average solution) to 1 (rarest solution) and is computed as twice 
	 * the absolute value of the difference between the CDF of the provided value in the adjusted 
	 * distribution and 0.5.
	 * @param value The value whose rarity is to be computed
	 * @return The rarity of the provided value
	 */
	public double computeRarity(double value)
	{
		if(distribution == null)
			return Double.NaN;
		
		if(value < min || value > max)
			return 1.0;
		
		return 2*Math.abs(0.5 - distribution.getCDF(value));
	}
	
}