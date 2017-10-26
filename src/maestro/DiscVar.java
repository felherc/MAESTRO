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

import probDist.DiscProbDist;
import utilities.Utilities;
import utilities.stat.ContSeries;

/**
 * Stores information for a discrete variable of the problem, including its name, the set of 
 * possible values and whether or not the variable corresponds to a scale.
 * @author Felipe Hernández
 */
public class DiscVar 
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the variable
	 */
	private String name;
	
	/**
	 * The minimum value the variable can take. This attribute and <code>count</code> are used when 
	 * the variable represents an actual number and not just an identifier. Conversely, the 
	 * attribute <code>values</code> is used when the variable represents identifiers. If the 
	 * attribute <code>values</code> is not null, the identifier representation is used.
	 */
	private int min;
	
	/**
	 * The number of possible values the variable can take. This attribute and <code>min</code> 
	 * are used when the variable represents an actual number and not just an identifier. 
	 * Conversely, the attribute <code>values</code> is used when the variable represents 
	 * identifiers. If the attribute <code>values</code> is not null, the identifier representation 
	 * is used.
	 */
	private int count;
	
	/**
	 * The identifiers of the values the variable can take. This attribute is used when the values 
	 * do not represent actual numbers but simply identifiers. Use the attributes <code>min</code> 
	 * and <code>count</code> to represent actual numbers instead. The instance is null when the 
	 * latter representation is used.
	 */
	private ArrayList<String> values;
	
	/**
	 * True if the values of the variable correspond to a scale
	 */
	private boolean scalar;
	
	/**
	 * The probability distribution that fits the values of the solutions in the population set. 
	 * The values are those corresponding to this variable. The distribution is used to compute the 
	 * rarity of new solutions.
	 */
	private DiscProbDist distribution;
	
	// --------------------------------------------------------------------------------------------
	// Constructors
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are the specified 
	 * number of integers starting from the <code>min</code> value and sequentially on. By default, 
	 * the variable is assumed to be scalar.
	 * @param name The identifier of the variable
	 * @param min The minimum integer value the variable can take
	 * @param count The number of values the variable can take
	 */
	public DiscVar(String name, int min, int count)
	{
		this.name = name;
		this.count = count;
		this.min = min;
		values = null;
		scalar = true;
		distribution = null;
	}
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are integers from 0 
	 * and on, that map to a set of identifiers. By default, the variable is assumed not to be 
	 * scalar.
	 * @param name The identifier of the variable
	 * @param values The set of identifiers of the variable's possible values
	 */
	public DiscVar(String name, ArrayList<String> values)
	{
		this.name = name;
		this.values = values;
		this.min = 0;
		this.count = values.size();
		scalar = false;
		distribution = null;
	}
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are the specified 
	 * number of integers starting from the <code>min</code> value and sequentially on.
	 * @param name The identifier of the variable
	 * @param min The minimum integer value the variable can take
	 * @param count The number of values the variable can take
	 * @param scalar True if the values of the variable correspond to a scale
	 */
	public DiscVar(String name, int min, int count, boolean scalar)
	{
		this.name = name;
		this.count = count;
		this.min = min;
		values = null;
		this.scalar = scalar;
		distribution = null;
	}
	
	/**
	 * Initializes a new discrete variable with its possible values. The values are integers from 0 
	 * and on, that map to the set of identifiers.
	 * @param name The identifier of the variable
	 * @param values The set of identifiers of the variable's possible values
	 * @param scalar True if the values of the variable correspond to a scale
	 */
	public DiscVar(String name, ArrayList<String> values, boolean scalar)
	{
		this.name = name;
		this.values = values;
		this.min = 0;
		this.count = values.size();
		this.scalar = scalar;
		distribution = null;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * @return The identifier of the variable
	 */
	public String getName() 
	{
		return name;
	}

	/**
	 * @return The minimum value the variable can take
	 */
	public int getMin() 
	{
		return min;
	}
	
	/**
	 * @return The maximum value the variable can take
	 */
	public int getMax() 
	{
		return min + count - 1;
	}

	/**
	 * @return The number of possible values the variable can take
	 */
	public int getCount() 
	{
		return count;
	}
	
	/**
	 * @return The value range of the variable (maximum value - minimum value)
	 */
	public double getRange()
	{
		return count - min;
	}

	/**
	 * @return The identifiers of the values the variable can take. Null if the values represent
	 * actual numbers.
	 */
	public ArrayList<String> getValues() 
	{
		return values;
	}

	/**
	 * @return True if the values of the variable correspond to a scale
	 */
	public boolean isScalar() 
	{
		return scalar;
	}

	/**
	 * @return The probability distribution that fits this variable values of the solutions in the 
	 * population set in order to compute the rarity of new solutions
	 */
	public DiscProbDist getDistribution() 
	{
		return distribution;
	}
	
	/**
	 * Generates random values to initialize the population
	 * @param count The number of values (solutions) to generate
	 * @return The random values generated
	 */
	public ArrayList<Integer> generateRandomValues(int count)
	{
		ArrayList<Integer> values = new ArrayList<Integer>();
		for(int i = 0 ; i < count ; i ++)
			values.add(min + Utilities.uniformRandomSelect(this.count));
		return values;
	}
	
	/**
	 * Validates that a candidate value is between {@link #min} and {@link #getMax}
	 * @param candidate The value to validate
	 * @return The validated value
	 */
	public int validate(int candidate) 
	{
		int validated	= candidate < min ? min : candidate;
		int max			= getMax();
		validated		= validated > max ? max : validated;
		return validated;
	}
	
	/**
	 * Adjust a probability distribution to the values observed in order to be able to compute the
	 * rarity of new solutions. If the variable is scalar, the adjusted distribution is continuous
	 * and the discrete values are obtained through rounded from it. If the variable is not scalar,
	 * a look-up table is used to store the probability of each value.
	 * @param values The values of the solutions in the current population set that correspond to 
	 * the variable
	 */
	public void adjustDistribution(ContSeries values)
	{
		if(scalar)
			distribution = new DiscProbDist(min, min + count - 1, values);
		else
			distribution = new DiscProbDist(values);
	}
	
	/**
	 * Computes the rarity of the provided value. The rarity measures how different or rare is the 
	 * value in comparison with the values in the first partition of the current population. The 
	 * rarity is a value from 0 (average solution) to 1 (rarest solution). If the variable is 
	 * scalar, the rarity is computed as twice the absolute value of the difference between the CDF 
	 * of the provided value in the adjusted distribution and 0.5. If the variable is not scalar, 
	 * the rarity is computed as the ratio between the probability of the value and the maximum 
	 * probability of all the values.
	 * @param value The value whose rarity is to be computed
	 * @return The rarity of the provided value
	 */
	public double computeRarity(int value)
	{
		if(distribution == null)
			return Double.NaN;
		
		if(value < min || value > min + count - 1)
			return 1.0;
		
		if(scalar)
			return 2*Math.abs(0.5 - distribution.getContinuous().getCDF(value));
		else
		{
			double maxProb = Collections.max(distribution.getProbs())/distribution.getProbSum();
			return 1 - distribution.getProb(value)/maxProb;
		}
	}
	
}
